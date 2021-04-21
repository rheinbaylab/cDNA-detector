import logging
try:
    from . import global_para
except ImportError:
    import global_para

def set_logger():
    logging.basicConfig(format= '%(asctime)s - %(levelname)s - %(message)s', level=logging.DEBUG)
    logger = logging.getLogger(__name__)
    handler = logging.FileHandler(global_para.out_log,'w')
    logger.addHandler(handler)
    handler.setLevel(logging.DEBUG)
    logFormatter = logging.Formatter(
        '%(asctime)s - %(levelname)s - %(message)s')
    handler.setFormatter(logFormatter)
    return logger