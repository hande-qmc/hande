def parse_yaml(stream, overrides={}):
    import yaml
    import sys
    from autocmake.interpolate import interpolate

    try:
        config = yaml.load(stream, yaml.SafeLoader)
    except yaml.YAMLError as exc:
        print(exc)
        sys.exit(-1)

    for k in config:
        if k in overrides:
            config[k] = overrides[k]

    config = interpolate(config, config)
    return config


def test_parse_yaml():
    text = """foo: bar
this: that
var: '%(foo)'
list:
  - a: '%(foo)'
  - b: '%(foo)'
  - c: '%(foo)'"""

    assert parse_yaml(text) == {'foo': 'bar', 'this': 'that', 'var': 'bar',
                                'list': [{'a': 'bar'}, {'b': 'bar'}, {'c': 'bar'}]}
