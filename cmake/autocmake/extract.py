def extract_list(config, section):
    from collections import Iterable
    l = []
    if 'modules' in config:
        for module in config['modules']:
            for k, v in module.items():
                for x in v:
                    if section in x:
                        if isinstance(x[section], Iterable) and not isinstance(x[section], str):
                            for y in x[section]:
                                l.append(y)
                        else:
                            l.append(x[section])
    return l


def to_d(l):
    """
    Converts list of dicts to dict.
    """
    _d = {}
    for x in l:
        for k, v in x.items():
            _d[k] = v
    return _d


def test_to_d():
    l = [{'a': 'b'}, {'c': 'd'}]
    d = {'a': 'b', 'c': 'd'}
    assert to_d(l) == d


def to_l(x):
    """
    Converts list of dicts to dict.
    """
    if isinstance(x, str):
        return [x]
    else:
        return x


def test_to_l():
    assert to_l('foo') == ['foo']
    assert to_l(['foo', 'bar']) == ['foo', 'bar']
