import yaml
import collections


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
