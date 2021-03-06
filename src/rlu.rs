#[derive(Debug, PartialEq)]
pub enum PivotPolicy {
    NoDiagonalElement,
    NoPivoting,
    PartialPivoting,
    ThresholdPivoting,
}

pub struct Options {
    pub pivot_policy: PivotPolicy,
    pub pivot_threshold: f64,
    pub drop_threshold: f64,
    pub col_fill_ratio: f64,
    pub fill_ratio: f64,
    pub expand_ratio: f64,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            pivot_policy: PivotPolicy::PartialPivoting,
            pivot_threshold: 1.0,
            drop_threshold: 0.0,  // do not drop
            col_fill_ratio: -1.0, // do not limit column fill ratio
            fill_ratio: 4.0,
            expand_ratio: 1.2,
        }
    }
}
