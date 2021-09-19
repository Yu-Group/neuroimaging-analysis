import time

def call_and_timing(fn, *args, _fn_name: str=None, _verbose: bool=True, **kwargs):
    start = time.time()
    result = fn(*args, **kwargs)
    duration = time.time() - start
    if _verbose:
        if _fn_name is None:
            print(f'Duration of call: {duration:.3}s')
        else:
            print(f'Duration of "{_fn_name}" call: {duration:.3}s')
    return result, duration
