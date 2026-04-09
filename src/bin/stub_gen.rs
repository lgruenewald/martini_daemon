use pyo3_stub_gen::Result;
use martini_daemon::stub_info;

// run with cargo run --bin stub_gen

fn main() -> Result<()> {
    // `stub_info` is a function defined by `define_stub_info_gatherer!` macro.
    let stub = stub_info()?;
    stub.generate()?;
    Ok(())
}