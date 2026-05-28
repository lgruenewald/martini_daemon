// WIP - this file is currently not used, martini_daemon/__formats/top_traj.py is used

#![allow(dead_code, unused)]
use flate2::Compression;
use flate2::Crc;
use flate2::read::ZlibDecoder;
use flate2::write::ZlibEncoder;
use pyo3::types::PyList;
use std::fs::File;
use std::fs::OpenOptions;
use std::io::{Read, Seek, SeekFrom, Write};
use std::u32;

use pyo3::{exceptions::PyException, prelude::*};

const MAGIC: [u8; 8] = [0xc0, b'T', b'O', b'P', b'T', b'R', 1, 0];

enum BufferState {
    Empty,
    HeaderWritten,
    AtomsWritten,
    BondsWritten,
}

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
    fn new(
        path: &str,
        title: &str,
        initial_molecules: Vec<(String, usize, usize)>,
        res_names: Vec<String>,
        res_ids: Vec<usize>,
        append: bool,
        truncate: usize,
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
            handle.seek(SeekFrom::Start(0))?;
            let mut magic = [0u8; 8];
            handle.read_exact(&mut magic)?;
            if magic != MAGIC {
                return Err(PyException::new_err(
                    "Magic number mismatch. Can't append to a different file format.",
                ));
            }

            // verify header
            // truncate if needed
            todo!()
        };
        Ok(Self {
            buffer: Vec::new(),
            buffer_state: BufferState::Empty,
            next_frame,
            n_atoms,
            handle: Some(handle),
        })
    }
}

#[pyclass]
pub struct TopTrajFrame {
    content: Vec<u8>,
    #[pyo3(get)]
    frame_index: u32,
    #[pyo3(get)]
    sim_step: u64,
    #[pyo3(get)]
    sim_time: f32,
    #[pyo3(get)]
    n_atoms: u32,

    // during initial decompression, only the non-py versions are filled
    // on first python access, python-owned obj's are constructed and the former is dropped
    // this should save a bit of memory and time, as often one or two fields only are ever read.
    // FIXME: more lazy access
    names_py: Option<Py<PyList>>,
    names: Option<Vec<String>>,
    types_py: Option<Py<PyList>>,
    types: Option<Vec<String>>,
    charges_py: Option<Py<PyList>>,
    charges: Option<Vec<f32>>,
    masses_py: Option<Py<PyList>>,
    masses: Option<Vec<f32>>,
    bonds_py: Option<Py<PyList>>,
    bonds: Option<Vec<(u32, u32)>>,
    // TODO continue here
}

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
}

#[pymethods]
impl TopTrajReader {
    #[new]
    fn new(path: String, py: Python<'_>) -> PyResult<Self> {
        let mut handle = std::fs::File::open(path.clone())?;

        handle.seek(SeekFrom::Start(0))?;
        let mut magic = [0u8; 8];
        handle.read_exact(&mut magic)?;
        if magic != MAGIC {
            return Err(PyException::new_err(
                "Magic number mismatch. Can't append to a different file format.",
            ));
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
        let start = chunk.start;
        let comp_len = chunk.comp_len;
        seek_to_chunk_end(start, comp_len, Some(expected_crc), &mut handle)?;

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
        })
    }
}
