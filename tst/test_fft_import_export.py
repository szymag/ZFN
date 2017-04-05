import logging
import os
import unittest

import numpy as np
<<<<<<< HEAD

from src.eig_problem.LoadFFT import LoadFFT
from src.fft_from_image.FFT import FFTfromImage
=======
from src.eig_problem.FFTfromFile import FFTfromFile

from src.fft_from_image.FFT import FFT
>>>>>>> 1D

logging.basicConfig(format='%(levelname)s:%(message)s',
                    filename='./log/test_fft_export_import.log',
                    filemode='w',
                    level=logging.DEBUG)


class TestFFTExportImport(unittest.TestCase):

    fft = FFTfromImage()
    fft_from_file = LoadFFT(9, 'I')

    def test_export_import(self):
        logging.info("### Running test_export_import...")
        lists_to_export = self.fft.wywolaj_fft2d(os.path.abspath("./tst/"))
        logging.info("numer of lists to export = %d" % len(lists_to_export))
        files = self.fft.wypisz_do_pliku(os.path.abspath("./tst/tmp/"), lists_to_export)
        logging.info("Lists were exported to files:")
        for filename in files:
            logging.info(" - %s" % filename)

        for i in range(len(files)):
            imported_list = self.fft_from_file.fourier_coefficient(files[i])
            np.testing.assert_array_almost_equal(imported_list, lists_to_export[i])
        logging.info("### end of test_export_import")

