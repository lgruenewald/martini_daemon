// Modified based on https://github.com/frobware/build-rs-example distributed under
// the MIT License.

// Copyright (c) 2025 Andrew McDermott

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.


//! `build.rs` - Generates a full version string at compile time.
//!
//! ## Purpose:
//! This script extracts version information (Git, Rust version, build date)
//! at compile time and sets it as `BUILD_VERSION`.
//!
//! ## Behaviour:
//! - Uses `CARGO_PKG_VERSION` as the base version.
//! - If Git is available, it appends:
//!   - The latest commit hash (`git rev-parse --short=10 HEAD`).
//!   - The dirty state if there are uncommitted changes.
//! - If Git is unavailable (e.g., tarball builds), Git info is omitted.
//! - The build date and Rust version are always included.

use std::{env, process::Command};

use chrono::Utc;

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    set_build_version();
}

/// Gets the Git commit hash and dirty status.
/// Returns `None` if Git is unavailable.
fn get_git_version() -> Option<String> {
    let commit_hash = Command::new("git")
        .args(["rev-parse", "--short=10", "HEAD"])
        .output()
        .ok()
        .filter(|output| output.status.success())
        .and_then(|output| {
            let commit_hash = String::from_utf8(output.stdout)
                .ok()
                .map(|s| s.trim().to_string());

            commit_hash
        });

    let is_dirty = Command::new("git")
        .args(["status", "--porcelain"])
        .output()
        .map(|out| !out.stdout.is_empty())
        .unwrap_or(false);

    let dirty_suffix = if is_dirty { "-dirty" } else { "" };
    Some(format!("{}{}", commit_hash?, dirty_suffix))
}

/// Sets `BUILD_VERSION`, the single authoritative version string.
fn set_build_version() {
    let build_date = Utc::now().format("%Y-%m-%dT%H:%M:%SZ").to_string();
    let git_version = get_git_version();
    let version = env!("CARGO_PKG_VERSION").to_string();

    let mut full_version = version.clone();

    full_version.push_str(" (");

    if let Some(git) = git_version {
        full_version.push_str(&format!("{git} "));
    }

    full_version.push_str(&format!("{build_date})"));
    
    println!("cargo:rustc-env=BUILD_VERSION={}", full_version);
}

