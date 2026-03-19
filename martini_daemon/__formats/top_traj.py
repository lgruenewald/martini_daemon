# Topology Trajectory format writer

class TopTrajWriter:
    """
    The Topology-Trajectory File Format .toptraj is described here.

    Goals of this format:
    - Store per-trajectory data for:
        - Simulation name
        - LIST OF INITIAL MOLECULES - a list of molecule names and counts. Only one per trajectory. Stored because
          it can be of interest for some analysis scripts.
    - Store per-frame data for:
        - Frame number, simulation step, simulation time, number of atoms
        - PER-ATOM DATA: atom names, residue names, residue IDs, atom types, charges, masses
        - LIST OF BONDS - all bonds, constraints, virtual sites. No data on the bond strengths, as the connectivity
          network is the most important, bond strength is a force field detail. Virtual sites always have bonds between
          the constructing atoms and the virtual particle, no extra bonds added among constructing atoms.
          It also sets a maximum of one bond between each pair of atoms, as the number of bonds as defined in a force
          field topology may not be so relevant for a connectivity graph.
        - LIST OF FRAGMENTS - list of dynamic atom groupings. a list of fragment name, fragment ID and constructing
          particles. In Martini Daemon, this stores the fragments found by the graph matching algorithm at each
          frame.
    - Store per-frame data as diffs from the previous frame
    - Optimize disk space and sequential access (both reading and writing)
    - Backward compatibility (it is desirable to be able to analyze data in old files), but no forward compatibility

    Assumptions made about the simulation:
        - less than 2^32 frames in total (the MD step can exceed 2^32 though)
        - less than 2^32 atoms
        - less than 2^32 residues
        - less than 2^32 fragments (by ID)

    Description of the binary layout:

    - TopTraj files consist of a header and body.
    - The header contains metadata only, and is present as raw bytes.
    - The body contains all the info described above. The header is compressed using ZLIB. ZLIB writes its own tiny
    header at the start of the body which contains the compression details. A compression level of 6 is recommended.
    - integers are all little endian and unsigned.
    - all strings are UTF-8 encoded and null terminated.
    - floats are 32 bit IEEE floats.
    - doubles are 64 bit IEEE floats.

    The header contains the following information:
    - The 6 byte Magic number "0xc0TOPTR"
    - 1 byte for major version of format (currently 1)
    - 1 byte for minor version of format (currently 0)

    The body contains the following information:
    - Simulation name
        - string - simulation title
    - List of initial molecules:
        - 4 byte integer - number of entries
        - List of entries:
            - string - molecule name
            - 4 byte integer - number of molecules in this entry
    - CRC32 checksum of header bytes+simulation name+list of initial molecules (as raw uncompressed bytes for the latter two)
    - List of frames:
        - 4 byte integer - frame number, must be one greater than previous frame, must be 0 for the first frame
        - 4 byte integer - number of atoms (n_atoms)
        - 8 byte integer - simulation step
        - double - simulation time, in nanoseconds
        - List of atom names
            - if first frame:
                - n_atoms strings
            - all other frames:
                - 4 byte integer n_changes (relative to previous frame)
                - for each change:
                    - 4 byte integer, index of changed atom
                    - string, new name
        - Similar setup for resname, resid, atom type, charges, mass. Payload type:
            - resname - string
            - resid - 4 byte integer
            - atom type - string
            - charge - float
            - mass - float
        - List of bonds:
            - in no particular order, however guaranteed to have at most one bond per i-j combination and no i-i
            - if first frame:
                - 8 byte integer n_bonds
                - for each bond:
                    - 4 byte integer, index of atom i
                    - 4 byte integer, index of atom j
            - all other frames:
                - 8 byte integer n_deletions (relative to previous frame)
                    - 4 byte integer, index of atom i
                    - 4 byte integer, index of atom j
                - readers should error if a bond is removed that was not there
                - 8 byte integer n_additions (relative to previous frame)
                - for each bond addition:
                    - 4 byte integer, index of atom i
                    - 4 byte integer, index of atom j
                - readers should error if a bond is added that was there
        - List of fragments:
            - if first frame:
                - 4 byte integer n_fragments
                - for each fragment:
                    - string, fragment name
                    - 4 byte integer, frag_id
                    - 4 byte integer, number of atoms
                        - n_frag_atoms 4 byte integers - containing atoms, if 0xFFFFFFFF it is a missing optional or forbidden atom
            - all other frames:
                - 4 byte integer n_deletions (relative to previous frame)
                    - 4 byte integer for FRAG_IDs of each deletion
                - 4 byte integer n_additions (relative to previous frame)
                    - string, fragment name
                    - 4 byte integer, frag_id
                    - 4 byte integer, number of atoms
                        - n_frag_atoms 4 byte integers - containing atoms, if 0xFFFFFFFF it is a missing optional or forbidden atom
        - CRC32 checksum of frame data
    """

    # TODO

