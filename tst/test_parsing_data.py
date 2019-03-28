import unittest
from src.io.DataReader import ParsingData


class ParsingDataTestCases(unittest.TestCase):
    parameter_set1 = ParsingData("./tst/Parameter_for_TheImpact.yaml")
    parameter_set2 = ParsingData("./tst/Sokolovskyy.yaml")

    @staticmethod
    def test_type():
        assert isinstance(ParsingDataTestCases.parameter_set1, ParsingData)
        assert isinstance(ParsingDataTestCases.parameter_set2, ParsingData)
        assert isinstance(ParsingDataTestCases.parameter_set1.input_parameters, dict)
        assert isinstance(ParsingDataTestCases.parameter_set2.input_parameters, dict)
        assert isinstance(ParsingData(ParsingDataTestCases.parameter_set1.input_parameters), ParsingData)
        assert isinstance(ParsingData(ParsingDataTestCases.parameter_set2.input_parameters), ParsingData)

    @staticmethod
    def test_parameter_type():
        assert ParsingDataTestCases.parameter_set2.rec_vector()[0] is None
        assert isinstance(ParsingDataTestCases.parameter_set1.physical_constant()[1], float)
        assert isinstance(ParsingDataTestCases.parameter_set1.fft_data(), str)

    @staticmethod
    def test_equality():
        assert ParsingDataTestCases.parameter_set1 != ParsingDataTestCases.parameter_set2
        assert ParsingDataTestCases.parameter_set2 == ParsingDataTestCases.parameter_set2

    @staticmethod
    def test_set_new_value():
        new_value1 = {'Mo': 1e6, 'l': 1e-17}
        new_value2 = {'mu0H0': 0.1, 'gamma': 194.6e9, 'mu0': 1e-6}
        new_value3 = 0.05

        ParsingDataTestCases.parameter_set1.set_new_value(new_value1, 'material_parameters', 'Co')
        assert ParsingDataTestCases.parameter_set1.material_constant('Co') == new_value1
        ParsingDataTestCases.parameter_set2.set_new_value(new_value2, 'physical_parameters')
        assert ParsingDataTestCases.parameter_set2.physical_constant() == (new_value2['gamma'], new_value2['mu0H0'])
        ParsingDataTestCases.parameter_set2.set_new_value(new_value3, 'physical_parameters', 'mu0H0')
        assert ParsingDataTestCases.parameter_set2.physical_constant()[1] == new_value3


if __name__ == '__main__':
    unittest.main()
