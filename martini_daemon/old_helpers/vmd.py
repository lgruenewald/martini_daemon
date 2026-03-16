# generate vmd bond info to visualize from .bonds and .xtc

# TODO move old_meta.py, gro_file.py and xtc file read/write to utils
# TODO give option for different xtc readers
import numpy as np
import zlib


xtc_readers = ["molly"]
try:
    import molly
    xtc_reader = "molly"
except ImportError:
    xtc_reader = None
    raise Warning(
        "Install one of the optional dependencies: "
        f"{', '.join(xtc_readers)} to generate vmd readable bonds."
    )


def generate_vmd_readable_bonds(
    frames: list[np.array], xtc_path: str, out_path: str
) -> None:
    """
    Given a list of bonds (from a .bonds file, loaded with
    old_reporters.bond_reporter.read_bonds), and a matching trajectory .xtc file
    (input) write a file in a format that is
    currently suitable for visualization by daemon.tcl.

    May take a while, reports progress to stdout.
    """

    if xtc_reader == "molly":
        xtc = molly.XTCReader(xtc_path)
    else:
        raise ImportError(
            "Generating vmd readable bonds requires one of the following "
            f"optional dependencies: {', '.join(xtc_readers)}."
        )
    print("Generating file for vmd bond visualization.")
    print(f"Using XTC reader {xtc_reader}.")

    comp_obj = zlib.compressobj(6)
    with open(out_path, "wb") as handle:
        for i, frame in enumerate(frames):
            if i % 8 == 0:
                print(f"Processing frame {i} of {len(frames)}")
            if xtc_reader == "molly":
                xtc_frame = xtc.pop_frame()
                pos = xtc_frame.positions
                box = np.array(
                    [
                        xtc_frame.box[0][0],
                        xtc_frame.box[1][1],
                        xtc_frame.box[2][2]
                    ],
                    dtype=np.float32
                )
            n_atoms = len(pos)
            bonds = [[] for _ in range(n_atoms)]
            crosses = np.any(
                np.abs(
                    pos[frame[:, 0]] - pos[frame[:, 1]]
                ) / box > 0.5,
                axis=-1
            )
            for [i, j], cross in zip(frame, crosses):
                if cross:
                    continue
                bonds[i].append(j)
                bonds[j].append(i)
            handle.write(
                comp_obj.compress(
                    (
                        "{" + "} {".join([
                            " ".join([
                                f"{other_atom}"
                                for other_atom in other_atom_list
                            ]) for other_atom_list in bonds
                        ]) + "} "
                    ).encode("utf-8")
                )
            )
        handle.write(comp_obj.flush())
    print("Generating file finished.")
