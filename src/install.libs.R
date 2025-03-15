# Check if cargo is installed
if (Sys.which("cargo") == "") {
  stop(
    "Rust (cargo) is required to build this package but is not installed.\n",
    "Please install Rust from https://www.rust-lang.org/ and try again."
  )
}
