#[cfg(feature = "debug")]
use std::io::Write;

#[cfg(feature = "debug")]
const MIN_WIDTH: usize = 5;

#[cfg(feature = "debug")]
const PADDING: usize = 1;

#[cfg(feature = "debug")]
const FLOAT_CONFIG: pretty_dtoa::FmtFloatConfig = pretty_dtoa::FmtFloatConfig::default()
    .add_point_zero(false)
    .max_significant_digits(6);

#[cfg(feature = "debug")]
pub fn matrix_table<I: crate::Int, S: crate::Scalar>(m: &crate::Matrix<I, S>) -> String {
    let mut tw = tabwriter::TabWriter::new(vec![])
        .minwidth(MIN_WIDTH)
        .padding(PADDING)
        .alignment(tabwriter::Alignment::Right);

    let n = m.len();

    for r in 0..n {
        for c in 0..n {
            match m[c].iter().find(|(i, _)| i.to_index() == r) {
                None => {
                    tw.write(b"-").unwrap();
                }
                Some((_, x)) => {
                    let s = x.pretty_string(FLOAT_CONFIG);
                    tw.write(s.as_bytes()).unwrap();
                }
            }
            if c == n - 1 {
                tw.write(b"\t\n").unwrap();
            } else {
                tw.write(b"\t").unwrap();
            }
        }
    }
    String::from_utf8(tw.into_inner().unwrap()).unwrap()
}

#[cfg(feature = "debug")]
macro_rules! debug {
    ($( $args:expr ),*) => { println!( $( $args ),* ); }
}

#[cfg(not(feature = "debug"))]
macro_rules! debug {
    ($( $args:expr ),*) => {};
}

pub(crate) use debug;
