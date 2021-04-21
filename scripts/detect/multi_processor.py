from joblib import Parallel,delayed
import multiprocessing
try:
    from . import global_para
except ImportError:
    import global_para

def apply_parallel(df_grouped, func):
    results = Parallel(n_jobs=global_para.n_threads,require='sharedmem')(delayed(func)(group) for name,group in df_grouped)
    return results

