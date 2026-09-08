// WIP - this file is currently not used, martini_daemon/__formats/top_traj.py is used

#![allow(dead_code, unused)]
use flate2::Compression;
use flate2::Crc;
use flate2::read::ZlibDecoder;
use flate2::write::ZlibEncoder;
use pyo3::exceptions::PyAssertionError;
use pyo3::exceptions::PyOSError;
use pyo3::exceptions::PyValueError;
use pyo3::ffi::PyExc_ZeroDivisionError;
use pyo3::types::PyFloat;
use pyo3::types::PyInt;
use pyo3::types::PyIterator;
use pyo3::types::PyString;
use pyo3::types::{PyAny, PyList, PySequence, PyType};
use pyo3_stub_gen::derive::gen_stub_pyfunction;
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};
use std::fs::File;
use std::fs::OpenOptions;
use std::io::{Read, Seek, SeekFrom, Write};
use std::sync::mpsc::RecvTimeoutError;
use std::u32;

use pyo3::{exceptions::PyException, prelude::*};

use crate::bond_graph::BondGraph;

const MAGIC: [u8; 8] = [0xc0, b'T', b'O', b'P', b'T', b'R', 1, 0];

#[derive(Debug, Copy, Clone)]
enum BufferState {
    Empty,
    HeaderWritten,
    AtomsWritten,
    BondsWritten,
}

#[gen_stub_pyclass]
#[pyclass]
pub struct TopTrajWriter {
    buffer: Vec<u8>,
    buffer_state: BufferState,
    next_frame: usize,
    n_atoms: usize,
    handle: Option<std::fs::File>,
}

impl<T> WriteExt for T
where
    T: Write,
{
    fn write_str(&mut self, value: &str) -> std::io::Result<usize> {
        let val_bytes = value.as_bytes();
        if val_bytes.len() <= 255 {
            self.write(&[val_bytes.len() as u8])?;
            self.write(val_bytes)
        } else {
            // number of bytes we can write
            let trunc = value.floor_char_boundary(255);
            assert!(trunc <= 255);
            self.write(&[trunc as u8])?;
            self.write(&val_bytes[0..trunc])
        }
    }

    fn write_u32(&mut self, value: u32) -> std::io::Result<usize> {
        self.write(&u32::to_le_bytes(value))
    }

    fn write_u64(&mut self, value: u64) -> std::io::Result<usize> {
        self.write(&u64::to_le_bytes(value))
    }

    fn write_f32(&mut self, value: f32) -> std::io::Result<usize> {
        self.write(&f32::to_le_bytes(value))
    }

    fn write_f64(&mut self, value: f64) -> std::io::Result<usize> {
        self.write(&f64::to_le_bytes(value))
    }

    fn write_chunk(&mut self, raw: &[u8]) -> std::io::Result<usize> {
        let mut crc = Crc::new();
        crc.update(raw);
        let amount = crc.amount();
        let crc = crc.sum();

        // FIXME: better buffering with less allocations
        // this mirrors the old python impl
        // API change idea: return a Chunk Writer
        let mut comp_buf = Vec::<u8>::new();
        let mut z = ZlibEncoder::new(&mut comp_buf, Compression::new(6));
        z.write(raw)?;
        z.try_finish()?;
        z.finish()?;
        self.write_u64(comp_buf.len() as u64)?;
        self.write_u64(raw.len() as u64)?;
        self.write(&comp_buf)?;
        assert!(raw.len() > u32::MAX as usize || amount == raw.len() as u32);
        self.write_u32(crc)
    }
}

trait WriteExt: Write {
    /// Write a string, with a single byte length prefix.
    fn write_str(&mut self, value: &str) -> std::io::Result<usize>;
    /// Write a little endian u32.
    fn write_u32(&mut self, value: u32) -> std::io::Result<usize>;
    /// Write a little endian u64.
    fn write_u64(&mut self, value: u64) -> std::io::Result<usize>;
    /// Write a little endian f32.
    fn write_f32(&mut self, value: f32) -> std::io::Result<usize>;
    /// Write a little endian f64.
    fn write_f64(&mut self, value: f64) -> std::io::Result<usize>;
    /// Write a compressed chunk.
    fn write_chunk(&mut self, raw: &[u8]) -> std::io::Result<usize>;
}

struct ChunkReader<'a> {
    inner: ZlibDecoder<&'a mut dyn Read>,
    start: u64,
    raw_len: u64,
    comp_len: u64,
    crc32: Crc,
}

impl<'a> ChunkReader<'a> {
    fn get_content(&mut self) -> std::io::Result<Vec<u8>> {
        let mut res: Vec<u8> = (0..self.raw_len).map(|_| 0u8).collect();

        self.inner.read_exact(&mut res)?;

        Ok(res)
    }

    fn get_crc(&mut self) -> std::io::Result<u32> {
        if self.inner.total_out() != self.raw_len || self.inner.total_in() != self.comp_len {
            return Err(std::io::Error::new(
                std::io::ErrorKind::Other,
                "Not all data in chunk was read. Chunk header possibly invalid.",
            ));
        }

        assert!(self.raw_len > u32::MAX as u64 || self.crc32.amount() == self.raw_len as u32);
        Ok(self.crc32.sum())
    }
}

/// Seek file to start + compressed len, read crc32 from file and compare with passed crc32, if provided.
///
/// Note: the seeking happens, so we don't make assumptions on how much the uncompressing algo read
fn seek_to_chunk_end(
    start: u64,
    comp_len: u64,
    crc32: Option<u32>,
    f: &mut std::fs::File,
) -> std::io::Result<()> {
    f.seek(SeekFrom::Start(start + comp_len + 8 + 8))?;
    let got_crc = f.read_u32()?;
    if crc32.is_some_and(|crc32| got_crc != crc32) {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "CRC32 Mismatch, the file is likely corrupt.",
        ));
    }
    Ok(())
}

impl<'a> Read for ChunkReader<'a> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let buf = if self.inner.total_out() + buf.len() as u64 > self.raw_len {
            &mut buf[0..(self.raw_len - self.inner.total_out()) as usize]
        } else {
            buf
        };

        // FIXME: .take() could be more elegant
        assert!(self.inner.total_out() <= self.raw_len);

        let n = self.inner.read(buf)?;

        if n > 0 {
            self.crc32.update(&buf[0..n]);
        }

        Ok(n)
    }
}

impl<T> ReadExt for T
where
    T: Read,
{
    fn read_str(&mut self) -> std::io::Result<String> {
        let mut len = [0u8; 1];
        self.read_exact(&mut len)?;
        let len = len[0];
        // FIXME: check if this does a useless allocation
        let mut res: Vec<u8> = (0..len).map(|_| 0u8).collect();
        self.read_exact(&mut res)?;

        Ok(String::from_utf8_lossy(&res).to_string())
    }

    fn read_u32(&mut self) -> std::io::Result<u32> {
        let mut buf = [0u8; 4];
        self.read_exact(&mut buf)?;
        Ok(u32::from_le_bytes(buf))
    }

    fn read_u64(&mut self) -> std::io::Result<u64> {
        let mut buf = [0u8; 8];
        self.read_exact(&mut buf)?;
        Ok(u64::from_le_bytes(buf))
    }

    fn read_f32(&mut self) -> std::io::Result<f32> {
        let mut buf = [0u8; 4];
        self.read_exact(&mut buf)?;
        Ok(f32::from_le_bytes(buf))
    }

    fn read_f64(&mut self) -> std::io::Result<f64> {
        let mut buf = [0u8; 8];
        self.read_exact(&mut buf)?;
        Ok(f64::from_le_bytes(buf))
    }
}

trait ReadExt: Read {
    fn read_str(&mut self) -> std::io::Result<String>;
    fn read_u32(&mut self) -> std::io::Result<u32>;
    fn read_u64(&mut self) -> std::io::Result<u64>;
    fn read_f32(&mut self) -> std::io::Result<f32>;
    fn read_f64(&mut self) -> std::io::Result<f64>;
}

/// Get a reader to read the inner chunk.
fn get_chunk_reader<'a>(reader: &'a mut std::fs::File) -> std::io::Result<ChunkReader<'a>> {
    let start = reader.seek(SeekFrom::Current(0))?;
    let comp_len = reader.read_u64()?;
    let raw_len = reader.read_u64()?;

    Ok(ChunkReader {
        inner: ZlibDecoder::new(reader),
        start,
        raw_len,
        comp_len,
        crc32: Crc::new(),
    })
}

// helpers
impl TopTrajWriter {
    fn write_frame_bonds_bond_list(&mut self, bonds: Vec<(u32, u32)>) -> PyResult<()> {
        // FIXME optimize
        let mut bond_graph = BondGraph::new(self.n_atoms);

        for (i, j) in bonds {
            bond_graph.add_bond(i as usize, j as usize);
        }
        self.write_frame_bonds_bond_graph(&bond_graph)
    }

    fn write_frame_bonds_bond_graph<'py>(&mut self, bonds: &BondGraph) -> PyResult<()> {
        match self.buffer_state {
            BufferState::AtomsWritten => (),
            _ => {
                return Err(PyException::new_err(
                    "write_frame_bonds() must be only called after calling write_frame_atoms().",
                ));
            }
        }
        let bonds = bonds.to_list();
        self.buffer.write_u64(bonds.len() as u64);

        for (i, j) in bonds {
            self.buffer.write_u32(i as u32);
            self.buffer.write_u32(j as u32);
        }
        self.buffer_state = BufferState::BondsWritten;

        Ok(())
    }
}

// python API
#[gen_stub_pymethods]
#[pymethods]
impl TopTrajWriter {
    /// Create a Topology Trajectory Writer.
    ///
    /// Note: this class owns a file handle. Call .finish() or use a context manager!
    ///
    /// :param path: The path to write to.
    /// :param title: The simulation title to write.
    /// :param initial_molecules: List of (name, count, atoms_per_mol) tuples.
    /// :param res_names: List of n_atoms residue names.
    /// :param res_ids: List of n_atoms residue IDs.
    /// :param append: Whether to append to an existing .toptraj file.
    /// :param truncate: If appending, the writer can optionally truncate the pre-existing
    ///   file. Specify the last MD step to keep.
    #[new]
    #[pyo3(signature=(path, title, initial_molecules, res_names, res_ids, /, append=false, truncate=None))]
    fn new(
        path: &str,
        title: &str,
        initial_molecules: Vec<(String, usize, usize)>,
        res_names: Vec<String>,
        res_ids: Vec<usize>,
        append: bool,
        truncate: Option<u64>,
    ) -> PyResult<Self> {
        let n_atoms = res_names.len();
        if res_names.len() != res_ids.len() {
            return Err(PyException::new_err("len(res_names) != len(res_ids)."));
        }

        let (handle, next_frame) = if !append {
            let mut handle = std::fs::File::create(path)?;
            let mut buf = Vec::<u8>::new();
            // write header
            handle.write(&MAGIC)?;

            buf.write_str(title)?;

            buf.write_u32(u32::try_from(initial_molecules.len())?)?;

            for (name, count, atoms_per_mol) in initial_molecules.iter() {
                buf.write_str(name)?;
                buf.write_u32(u32::try_from(*count)?)?;
                buf.write_u32(u32::try_from(*atoms_per_mol)?)?;
            }

            buf.write_u32(u32::try_from(n_atoms)?)?;

            for name in res_names.iter() {
                buf.write_str(name)?;
            }

            for id in res_ids.iter() {
                buf.write_u32(u32::try_from(*id)?)?;
            }

            handle.write_chunk(&buf)?;
            handle.flush()?;

            (handle, 0)
        } else {
            let mut handle = OpenOptions::new().read(true).write(true).open(path)?;
            let file_end = handle.seek(SeekFrom::End(0))?;
            handle.seek(SeekFrom::Start(0))?;
            let mut magic = [0u8; 8];
            handle.read_exact(&mut magic)?;
            if magic != MAGIC {
                return Err(PyException::new_err(
                    "Magic number mismatch. Can't append to a different file format.",
                ));
            }

            // truncate if needed
            let chunk = get_chunk_reader(&mut handle)?;
            let mut n_frames = 0;
            // FIXME do more header verifications
            seek_to_chunk_end(chunk.start, chunk.comp_len, None, &mut handle);

            loop {
                let start = handle.seek(SeekFrom::Current(0))?;

                if start >= file_end {
                    // EOF -> break
                    break;
                }

                // decompress chunk header
                let mut chunk = get_chunk_reader(&mut handle)?;
                let frame_index = chunk.read_u32()?;
                let n_atoms = chunk.read_u32()?;
                let sim_step = chunk.read_u64()?;
                let sim_time = chunk.read_f64()?;

                if let Some(truncate) = truncate
                    && truncate > sim_step
                {
                    // ignore frame, go back to before
                    handle.seek(SeekFrom::Start(start));
                    break;
                } else {
                    // accept frame, go to the end
                    n_frames += 1;
                    seek_to_chunk_end(chunk.start, chunk.comp_len, None, &mut handle);
                }
            }

            let curr = handle.seek(SeekFrom::Current(0))?;
            handle.set_len(curr)?;

            (handle, n_frames)
        };
        Ok(Self {
            buffer: Vec::new(),
            buffer_state: BufferState::Empty,
            next_frame,
            n_atoms,
            handle: Some(handle),
        })
    }

    fn __enter__(self_: PyRef<'_, Self>) -> PyRef<'_, Self> {
        self_
    }

    fn __exit__(
        &mut self,
        _exc_type: Option<&Bound<'_, PyType>>,
        _exc_val: Option<&Bound<'_, PyAny>>,
        _exc_tb: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<bool> {
        self.finish()?;
        // do not suppress exceptions
        Ok(false)
    }

    /// Close the writer.
    ///
    /// Identical to finish().
    fn close(&mut self) -> PyResult<()> {
        self.finish()
    }

    /// Close the writer.
    fn finish(&mut self) -> PyResult<()> {
        self.handle = None;
        Ok(())
    }

    // FIXME: in the future a Frame Builder would be nicer, but I don't want any breaking changes atm

    /// Write the frame header to the writer buffer.
    ///
    /// Frames should be written by subsequent calls to the Writer in this
    /// order:
    ///
    /// * new_frame()
    /// * write_frame_atoms()
    /// * write_frame_bonds()
    /// * write_frame()
    ///
    /// :param frame_num: The current frame number. Must be one larger than the previous frame.
    /// :param sim_step: The simulation step.
    /// :param time_ps: The simulation time in picoseconds.
    /// :param n_atoms: The number of atoms in this frame.
    fn new_frame(
        &mut self,
        frame_num: u32,
        sim_step: u64,
        time_ps: f64,
        n_atoms: u32,
    ) -> PyResult<()> {
        match self.buffer_state {
            BufferState::Empty => (),
            _ => {
                return Err(PyException::new_err(
                    "new_frame() must be only called after creating a new writer, or after write_frame().",
                ));
            }
        }

        if self.next_frame != frame_num as usize {
            return Err(PyValueError::new_err(
                "Specified frame_num must be 1 larger than the previous frame.",
            ));
        }

        if self.n_atoms != n_atoms as usize {
            return Err(PyValueError::new_err(
                "n_atoms in frame and on Writer construction not identical.",
            ));
        }

        assert!(self.buffer.len() == 0);

        self.buffer.write_u32(frame_num)?;
        self.buffer.write_u32(n_atoms)?;
        self.buffer.write_u64(sim_step)?;
        self.buffer.write_f64(time_ps)?;
        self.buffer_state = BufferState::HeaderWritten;

        Ok(())
    }

    /// Write the current frame atom information to disk.
    ///
    /// The number of atoms must be the same as n_atoms specified in new_frame().
    /// new_frame() must be called first. register_frame_atoms() must be called exactly once per frame.
    ///
    /// :param names: The names of the atoms.
    /// :param atom_types: The Non-Bonded force atom types, per atom.
    /// :param charges: The charges of the atoms, per atom.
    /// :param masses: The masses of the atoms, per atom.
    fn write_frame_atoms<'py>(
        &mut self,
        #[gen_stub(override_type(type_repr = "list[str]"))] names: Bound<'py, PySequence>,
        #[gen_stub(override_type(type_repr = "list[str]"))] atom_types: Bound<'py, PySequence>,
        #[gen_stub(override_type(type_repr = "list[float]"))] charges: Bound<'py, PySequence>,
        #[gen_stub(override_type(type_repr = "list[float]"))] masses: Bound<'py, PySequence>,
    ) -> PyResult<()> {
        match self.buffer_state {
            BufferState::HeaderWritten => (),
            _ => {
                return Err(PyException::new_err(
                    "write_frame_atoms() must be only called after calling new_frame().",
                ));
            }
        }
        let n_atoms = self.n_atoms;
        if names.len()? != n_atoms
            || atom_types.len()? != n_atoms
            || charges.len()? != n_atoms
            || masses.len()? != n_atoms
        {
            return Err(PyValueError::new_err(
                "All arguments to write_frame_atoms must have length n_atoms",
            ));
        }

        for i in 0..n_atoms {
            let Ok(name) = names.get_item(i)?.cast_into::<PyString>() else {
                return Err(PyValueError::new_err("All names must be of type str."));
            };
            let name: String = name.extract()?;
            self.buffer.write_str(&name);
        }
        for i in 0..n_atoms {
            let Ok(atom_type) = atom_types.get_item(i)?.cast_into::<PyString>() else {
                return Err(PyValueError::new_err("All atom types must be of type str."));
            };
            let atom_type: String = atom_type.extract()?;
            self.buffer.write_str(&atom_type);
        }
        for i in 0..n_atoms {
            let Ok(charge) = charges.get_item(i)?.cast_into::<PyFloat>() else {
                return Err(PyValueError::new_err("All charges must be of type float."));
            };
            let charge: f32 = charge.extract()?;
            self.buffer.write_f32(charge);
        }
        for i in 0..n_atoms {
            let Ok(mass) = masses.get_item(i)?.cast_into::<PyFloat>() else {
                return Err(PyValueError::new_err("All masses must be of type float."));
            };
            let mass: f32 = mass.extract()?;
            self.buffer.write_f32(mass);
        }
        self.buffer_state = BufferState::AtomsWritten;
        Ok(())
    }

    /// Write the bonds for the current frame to disk.
    ///
    /// Note: will automatically remove duplicates, self-bonds
    ///       and will re-order bonds to have smaller first in each entry.
    ///
    /// :param bonds: The bonds to write, as a BondGraph object or as list of (i, j) tuples.
    fn write_frame_bonds<'py>(&mut self, bonds: Bound<'py, PyAny>) -> PyResult<()> {
        if let Ok(bond_list) = bonds.extract::<Vec<(u32, u32)>>() {
            self.write_frame_bonds(bonds)
        } else if let Ok(bond_graph) = bonds.extract::<Bound<'py, BondGraph>>() {
            self.write_frame_bonds_bond_graph(&*bond_graph.borrow())
        } else {
            Err(PyValueError::new_err(
                "Bond list must be provided as a list of (i, j) tuples or a BondGraph object.",
            ))
        }
    }

    /// Finish writing the current frame to disk.
    ///
    /// Must call register_frame_atoms and register_frame_bonds exactly once first.
    ///
    /// Note: will automatically flush after finishing the frame.
    fn write_frame(&mut self) -> PyResult<()> {
        match self.buffer_state {
            BufferState::BondsWritten => (),
            _ => {
                return Err(PyException::new_err(
                    "write_frame() must be only called after calling write_frame_bonds().",
                ));
            }
        }

        let Some(handle) = &mut self.handle else {
            return Err(PyOSError::new_err("Writer is already closed."));
        };

        handle.write_chunk(&self.buffer);
        handle.flush();

        self.buffer.clear();
        self.buffer_state = BufferState::Empty;
        self.next_frame += 1;
        Ok(())
    }

    fn debug_print(&self) {
        println!("TopTrajWriter debugprint");
        println!("buffer: {:?}", self.buffer);
        println!("buffer state: {:?}", self.buffer_state);
        println!(
            "next_frame: {:?} n_atoms: {:?} handle: {:?}",
            self.next_frame, self.n_atoms, self.handle
        );
    }
}

#[gen_stub_pyclass]
#[pyclass]
pub struct TopTrajFrame {
    #[pyo3(get)]
    frame_index: u32,
    #[pyo3(get)]
    sim_step: u64,
    #[pyo3(get)]
    sim_time: f64,
    #[pyo3(get)]
    n_atoms: u32,

    // during initial decompression, only the non-py versions are filled
    // on first python access, python-owned obj's are constructed and the former is dropped
    // this should save a bit of memory and time, as often one or two fields only are ever read.
    // FIXME: more lazy access, numpy access, bond graph access
    names_py: Option<Py<PyList>>,
    names: Option<Vec<String>>,
    atom_types_py: Option<Py<PyList>>,
    atom_types: Option<Vec<String>>,
    charges_py: Option<Py<PyList>>,
    charges: Option<Vec<f32>>,
    masses_py: Option<Py<PyList>>,
    masses: Option<Vec<f32>>,
    bonds_py: Option<Py<PyList>>,
    bonds: Option<Vec<(u32, u32)>>,
}

#[gen_stub_pymethods]
#[pymethods]
impl TopTrajFrame {
    #[getter]
    fn names(&mut self, py: Python<'_>) -> PyResult<Py<PyList>> {
        match &mut self.names_py {
            Some(list) => {
                assert!(self.names == None);
                Ok(list.clone_ref(py))
            }
            None => {
                if let Some(names) = std::mem::take(&mut self.names) {
                    assert!(self.names == None);
                    self.names_py = Some(PyList::new(py, names)?.unbind());
                    Ok(self.names_py.as_mut().unwrap().clone_ref(py))
                } else {
                    Err(PyAssertionError::new_err(
                        "Internal error: no names in frame.",
                    ))
                }
            }
        }
    }
    #[getter]
    fn atom_types(&mut self, py: Python<'_>) -> PyResult<Py<PyList>> {
        match &mut self.atom_types_py {
            Some(list) => {
                assert!(self.atom_types == None);
                Ok(list.clone_ref(py))
            }
            None => {
                if let Some(types) = std::mem::take(&mut self.atom_types) {
                    assert!(self.atom_types == None);
                    self.atom_types_py = Some(PyList::new(py, types)?.unbind());
                    Ok(self.atom_types_py.as_mut().unwrap().clone_ref(py))
                } else {
                    Err(PyAssertionError::new_err(
                        "Internal error: no types in frame.",
                    ))
                }
            }
        }
    }
    #[getter]
    fn charges(&mut self, py: Python<'_>) -> PyResult<Py<PyList>> {
        match &mut self.charges_py {
            Some(list) => {
                assert!(self.charges == None);
                Ok(list.clone_ref(py))
            }
            None => {
                if let Some(charges) = std::mem::take(&mut self.charges) {
                    assert!(self.charges == None);
                    self.charges_py = Some(PyList::new(py, charges)?.unbind());
                    Ok(self.charges_py.as_mut().unwrap().clone_ref(py))
                } else {
                    Err(PyAssertionError::new_err(
                        "Internal error: no charges in frame.",
                    ))
                }
            }
        }
    }
    #[getter]
    fn masses(&mut self, py: Python<'_>) -> PyResult<Py<PyList>> {
        match &mut self.masses_py {
            Some(list) => {
                assert!(self.masses == None);
                Ok(list.clone_ref(py))
            }
            None => {
                if let Some(masses) = std::mem::take(&mut self.masses) {
                    assert!(self.masses == None);
                    self.masses_py = Some(PyList::new(py, masses)?.unbind());
                    Ok(self.masses_py.as_mut().unwrap().clone_ref(py))
                } else {
                    Err(PyAssertionError::new_err(
                        "Internal error: no masses in frame.",
                    ))
                }
            }
        }
    }
    #[getter]
    fn bonds(&mut self, py: Python<'_>) -> PyResult<Py<PyList>> {
        match &mut self.bonds_py {
            Some(list) => {
                assert!(self.bonds == None);
                Ok(list.clone_ref(py))
            }
            None => {
                if let Some(bonds) = std::mem::take(&mut self.bonds) {
                    assert!(self.bonds == None);
                    self.bonds_py = Some(PyList::new(py, bonds)?.unbind());
                    Ok(self.bonds_py.as_mut().unwrap().clone_ref(py))
                } else {
                    Err(PyAssertionError::new_err(
                        "Internal error: no bonds in frame.",
                    ))
                }
            }
        }
    }
}

#[gen_stub_pyclass]
#[pyclass]
pub struct TopTrajReader {
    #[pyo3(get)]
    path: String,
    #[pyo3(get)]
    major_version: u8,
    #[pyo3(get)]
    minor_version: u8,
    #[pyo3(get)]
    header: [u8; 8],
    #[pyo3(get)]
    title: String,
    #[pyo3(get)]
    initial_molecules: Py<PyList>,
    #[pyo3(get)]
    n_atoms: u32,
    #[pyo3(get)]
    res_names: Py<PyList>,
    #[pyo3(get)]
    res_ids: Py<PyList>,

    /// Offsets of the starts of all frames + the end of the file
    #[pyo3(get)]
    frame_offsets: Vec<u64>,

    handle: Option<std::fs::File>,
}

impl TopTrajReader {
    fn handle(&mut self) -> PyResult<&mut File> {
        match &mut self.handle {
            Some(handle) => Ok(handle),
            None => Err(PyOSError::new_err("Reader is closed.")),
        }
    }
}

#[gen_stub_pymethods]
#[pymethods]
impl TopTrajReader {
    #[new]
    /// Create a new TopTrajReader.
    fn new(path: String, py: Python<'_>) -> PyResult<Self> {
        let mut handle = std::fs::File::open(path.clone())?;

        let file_end = handle.seek(SeekFrom::End(0))?;

        handle.seek(SeekFrom::Start(0))?;
        let mut magic = [0u8; 8];
        handle.read_exact(&mut magic)?;
        if magic != MAGIC {
            return Err(PyException::new_err("Magic number mismatch."));
        }

        let mut chunk = get_chunk_reader(&mut handle)?;

        let title = chunk.read_str()?;

        // FIXME: pass iterator to PyList::new()
        let n_initial_molecules = chunk.read_u32()?;

        if n_initial_molecules > u16::MAX as u32 {
            // FIXME: this fixes allocation panics on bad files in a not very elegant manner
            return Err(PyException::new_err(
                "Initial molecules > 2**16, the file is likely corrupt.",
            ));
        }

        let mut initial_molecules: Vec<(String, u32, u32)> =
            Vec::with_capacity(n_initial_molecules as usize);
        for _ in 0..n_initial_molecules {
            initial_molecules.push((chunk.read_str()?, chunk.read_u32()?, chunk.read_u32()?));
        }

        let n_atoms: u32 = chunk.read_u32()?;

        let mut res_names: Vec<String> = Vec::with_capacity(n_atoms as usize);
        let mut res_ids: Vec<u32> = Vec::with_capacity(n_atoms as usize);

        for _ in 0..n_atoms {
            res_names.push(chunk.read_str()?);
        }

        for _ in 0..n_atoms {
            res_ids.push(chunk.read_u32()?);
        }

        let expected_crc = chunk.get_crc()?;
        seek_to_chunk_end(chunk.start, chunk.comp_len, Some(expected_crc), &mut handle)?;

        let mut offsets: Vec<u64> = Vec::new();

        // precalculate offsets for fast reading in the future
        loop {
            let start = handle.seek(SeekFrom::Current(0))?;

            if start >= file_end {
                break;
            }

            let mut chunk = get_chunk_reader(&mut handle)?;
            seek_to_chunk_end(chunk.start, chunk.comp_len, None, &mut handle)?;
            offsets.push(start);
        }

        assert!(handle.seek(SeekFrom::Current(0))? == file_end);
        handle.seek(SeekFrom::Start(offsets[0]));

        offsets.push(file_end);

        Ok(Self {
            path,
            major_version: 1,
            minor_version: 0,
            title,
            initial_molecules: PyList::new(py, initial_molecules)?.unbind(),
            header: magic,
            n_atoms,
            res_names: PyList::new(py, res_names)?.unbind(),
            res_ids: PyList::new(py, res_ids)?.unbind(),
            frame_offsets: offsets,
            handle: Some(handle),
        })
    }

    /// Return the current position in the file.
    ///
    /// :param frame_index: If specified, return the position of this frame in the file.
    #[pyo3(signature=(frame_index=None))]
    fn tell(&mut self, frame_index: Option<usize>) -> PyResult<u64> {
        match frame_index {
            None => Ok(self.handle()?.seek(SeekFrom::Current(0))?),
            Some(index) => {
                if self.frame_offsets.len() <= index {
                    return Err(PyValueError::new_err("Index out of range."));
                }
                Ok(self.frame_offsets[index])
            }
        }
    }

    // FIXME: change with seek_to(frame_index) instead, once I want breaking changes
    /// Set the reader to position, as returned by tell().
    fn seek(&mut self, pos: u64) -> PyResult<()> {
        self.handle()?.seek(SeekFrom::Start(pos))?;
        Ok(())
    }

    /// Read a frame from the toptraj file.
    ///
    /// :return: A toptraj frame, or None if finished.
    fn read_frame(&mut self) -> PyResult<Option<TopTrajFrame>> {
        if self.handle()?.seek(SeekFrom::Current(0))? >= *self.frame_offsets.last().unwrap() {
            return Ok(None);
        }
        let mut chunk = get_chunk_reader(&mut *self.handle()?)?;

        let frame_index = chunk.read_u32()?;
        let n_atoms = chunk.read_u32()?;
        let sim_step = chunk.read_u64()?;
        let sim_time = chunk.read_f64()?;

        let names = (0..n_atoms)
            .map(|_| chunk.read_str())
            .collect::<Result<Vec<String>, std::io::Error>>()?;
        let types = (0..n_atoms)
            .map(|_| chunk.read_str())
            .collect::<Result<Vec<String>, std::io::Error>>()?;
        let charges = (0..n_atoms)
            .map(|_| chunk.read_f32())
            .collect::<Result<Vec<f32>, std::io::Error>>()?;
        let masses = (0..n_atoms)
            .map(|_| chunk.read_f32())
            .collect::<Result<Vec<f32>, std::io::Error>>()?;
        let n_bonds = chunk.read_u64()?;
        let mut bonds = Vec::with_capacity(n_bonds as usize);

        for i in 0..n_bonds {
            bonds.push((chunk.read_u32()?, chunk.read_u32()?));
        }

        let mut frame = TopTrajFrame {
            frame_index,
            sim_step,
            sim_time,
            n_atoms,
            names: Some(names),
            names_py: None,
            atom_types: Some(types),
            atom_types_py: None,
            charges: Some(charges),
            charges_py: None,
            masses: Some(masses),
            masses_py: None,
            bonds: Some(bonds),
            bonds_py: None,
        };

        seek_to_chunk_end(
            chunk.start,
            chunk.comp_len,
            None,
            self.handle.as_mut().unwrap(),
        )?;

        Ok(Some(frame))
    }

    /// Skip a frame.
    ///
    /// If a frame was skipped, return True.
    /// If the reader is already at EOF, return False.
    fn skip_frame(&mut self) -> PyResult<bool> {
        if self.handle()?.seek(SeekFrom::Current(0))? >= *self.frame_offsets.last().unwrap() {
            return Ok(false);
        }
        let chunk = get_chunk_reader(&mut *self.handle()?)?;
        seek_to_chunk_end(
            chunk.start,
            chunk.comp_len,
            None,
            self.handle.as_mut().unwrap(),
        )?;
        Ok(true)
    }

    /// Close the reader.
    fn finish(&mut self) -> PyResult<()> {
        self.handle = None;
        Ok(())
    }

    /// Close the reader.
    ///
    /// Identical to finish().
    fn close(&mut self) -> PyResult<()> {
        self.finish()
    }

    fn __enter__(self_: PyRef<'_, Self>) -> PyRef<'_, Self> {
        self_
    }

    fn __exit__(
        &mut self,
        _exc_type: Option<&Bound<'_, PyType>>,
        _exc_val: Option<&Bound<'_, PyAny>>,
        _exc_tb: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<bool> {
        self.finish()?;
        // do not suppress exceptions
        Ok(false)
    }
}
