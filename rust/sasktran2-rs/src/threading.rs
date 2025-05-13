use std::sync::{Arc, Mutex};

use anyhow::Result;

static THREADPOOL: Mutex<Option<Arc<rayon::ThreadPool>>> = Mutex::new(None);

fn create_pool(num_threads: usize) -> Result<rayon::ThreadPool> {
    match rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
    {
        Err(e) => Err(e.into()),
        Ok(pool) => Ok(pool),
    }
}

pub fn set_num_threads(num_threads: usize) -> Result<()> {
    let mut pool = THREADPOOL.lock().unwrap();
    *pool = Some(create_pool(num_threads)?.into());
    Ok(())
}

pub fn thread_pool() -> Result<Arc<rayon::ThreadPool>, rayon::ThreadPoolBuildError> {
    let mut pool = THREADPOOL.lock().unwrap();

    if pool.is_none() {
        *pool = Some(Arc::new(
            rayon::ThreadPoolBuilder::new().num_threads(1).build()?,
        ));
    }

    Ok(pool.as_ref().unwrap().clone())
}
