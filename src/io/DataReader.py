import yaml
import warnings


class ParsingData:
    # TODO: Rename class, I guess
    # TODO: consider static method
    def __init__(self, input_parameters):
        if isinstance(input_parameters, str):
            self.input_parameters = self.load_yaml_file(input_parameters)
        elif isinstance(input_parameters, dict):
            self.input_parameters = input_parameters
        else:
            raise IOError('Input parameters not understood: ' + str(input_parameters))

    def __eq__(self, other):
        return self.input_parameters == other.input_parameters

    def _ne__(self, other):
        return self.input_parameters != other.input_parameters

    def load_yaml_file(self, input_parameters):
        with open(input_parameters, 'r') as stream:
            try:
                return self.convert_to_value(yaml.load(stream))
            except yaml.YAMLError as exc:
                return exc

    def convert_to_value(self, d):
        if len(d) == d or type(d) != dict:
            return d
        for k, v in d.items():
            if isinstance(v, dict):
                self.convert_to_value(v)
            elif isinstance(v, str):
                try:
                    d[k] = float(v)
                except Exception:
                    if d[k] == 'None':
                        d[k] = None
                    elif d[k] == 'True':
                        d[k] = True
                    elif d[k] == 'False':
                        d[k] = False
                    else:
                        pass
        return d

    def set_new_value(self, value, argument_1, argument_2=None):
        if argument_2 is None:
            if not isinstance(self.input_parameters[argument_1], type(value)):
                warnings.warn('Old values has ' + str(type(self.input_parameters[argument_1]))
                          + ' while new one has ' + str(type(value)))
            self.input_parameters[argument_1] = value
        else:

            if not isinstance(self.input_parameters[argument_1][argument_2], type(value)):
                warnings.warn('Old values has ' + str(type(self.input_parameters[argument_1][argument_2]))
                          + ' while new one has ' + str(type(value)))
            self.input_parameters[argument_1][argument_2] = value

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

    def output_file(self, data_type):
        if data_type == 'dispersion':
            return self.input_parameters['numerical_parameters']['output_file'] + '.dys'
        elif data_type == 'vectors':
            return self.input_parameters['numerical_parameters']['output_file'] + '.vec'
        else:
            return self.input_parameters['numerical_parameters']['output_file']

    def x(self):
        return self.input_parameters['system_dimensions']['x']

    def material_constant(self, material_name):
        return self.input_parameters['material_parameters'][material_name]

    def mu0(self):
        return self.input_parameters['physical_parameters']['mu0']

    def angle(self):
        return self.input_parameters['system_dimensions']['angle']

    def perpendicular_bloch_vector(self):
        return self.input_parameters['perp_vector']['angle'], self.input_parameters['perp_vector']['max_value']


if __name__ == "__main__":
    pass
