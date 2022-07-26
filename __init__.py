import logging
default_logger_format = ('%(levelname)s: [%(asctime)s] %(name)s'
                         ' - %(message)s')
default_date_format = '%Y-%m-%d %H:%M:%S'
# ms_logger_format = ('%(levelname)s: [%(asctime)s.%(msecs)03d] %(name)s'
#                          ' - %(message)s')
logging.basicConfig(format=default_logger_format, level=logging.INFO,#format=ms_logger_format,
                    datefmt=default_date_format)


# __version__ = '1.0.0'
