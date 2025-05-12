use anyhow::Result;

pub fn create_pool(num_threads: usize) -> Result<rayon::ThreadPool> {
    match rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
    {
        Err(e) => Err(e.into()),
        Ok(pool) => Ok(pool),
    }
}
