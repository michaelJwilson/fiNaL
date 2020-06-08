import yaml

def get_yaml(path):
  with open(path, 'r') as stream:
    result = yaml.safe_load(stream)

  return  result
