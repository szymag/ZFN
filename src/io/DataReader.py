import yaml


def to_float(d):
    if len(d) == d or type(d) != dict:
        return d
    for k, v in d.items():
        if isinstance(v, dict):
            to_float(v)
        elif isinstance(v, str):
            try:
                d[k] = float(v)
            except Exception:
                pass
    return d


def load_yaml_file(file_name):
    with open(file_name, 'r') as stream:
        try:
            return to_float(yaml.load(stream))
        except yaml.YAMLError as exc:
            return to_float(exc)


class ParsingData:
    # TODO: Rename class, I guess
    # TODO: put in above functions into this class
    # TODO: consider static method
    def __init__(self, input_parameters):
        if isinstance(input_parameters, str):
            self.input_parameters = load_yaml_file(input_parameters)
        elif isinstance(input_parameters, dict):
            # TODO: when data are loaded in GUI, they should be stored in dict
            self.input_parameters = input_parameters
        else:
            raise IOError

    def bloch_vector(self):
        return self.input_parameters['q_vector']['start'], self.input_parameters['q_vector']['end'],\
               self.input_parameters['q_vector']['dispersion_count']

    def lattice_const(self):
        return self.input_parameters['system_dimensions']['a'], self.input_parameters['system_dimensions']['b']

    def thickness(self):
        return self.input_parameters['system_dimensions']['d']

    def rec_vector(self):
        return self.input_parameters['numerical_parameters']['rec_vector_x'],\
               self.input_parameters['numerical_parameters']['rec_vector_y']

    def physical_constant(self):
        return self.input_parameters['physical_parameters']['gamma'],\
               self.input_parameters['physical_parameters']['mu0H0']

    def input_fft_file(self):
        return self.input_parameters['numerical_parameters']['fft_file']

    def output_file(self):
        if self.input_parameters['numerical_parameters']['output_file'][-4:] == '.vec':
            return self.input_parameters['numerical_parameters']['output_file']
        else:
            return self.input_parameters['numerical_parameters']['output_file'] + '.vec'

    def x(self):
        return self.input_parameters['system_dimensions']['x']

    def material_constant(self, material_name):
        return self.input_parameters['material_parameters'][material_name]

    def mu0(self):
        return self.input_parameters['physical_parameters']['mu0']

    def angle(self):
        return self.input_parameters['system_dimensions']['angle']