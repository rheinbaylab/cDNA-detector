import logging

def set_logger(out_log):
    logging.basicConfig(format= '%(asctime)s - %(levelname)s - %(message)s', level=logging.DEBUG)
    # logger = logging.getLogger(__name__).addHandler(logging.FileHandler('myLogs.log'))
    logging.logger = logging.getLogger(__name__)
    handler = logging.FileHandler(out_log,'w')
    logging.logger.addHandler(handler)
    handler.setLevel(logging.DEBUG)
    logFormatter = logging.Formatter(
        '%(asctime)s - %(levelname)s - %(message)s')
    handler.setFormatter(logFormatter)
