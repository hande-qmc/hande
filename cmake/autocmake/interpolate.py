def replace(s, d):
    from re import findall

    if isinstance(s, str):
        for var in findall(r"%\(([A-Za-z0-9_]*)\)", s):
            s = s.replace("%({0})".format(var), str(d[var]))
    return s


def test_replace():
    assert replace('hey %(foo) ho %(bar)',
                   {'foo': 'hey', 'bar': 'ho'}) == 'hey hey ho ho'


def interpolate(d, d_map):
    from collections import Mapping, Iterable
    from copy import copy

    for k, v in d.items():
        if isinstance(v, Mapping):
            d[k] = interpolate(d[k], d_map)
        elif isinstance(v, Iterable) and not isinstance(v, str):
            l = []
            for x in v:
                if isinstance(x, Mapping):
                    l.append(interpolate(x, d_map))
                else:
                    l.append(replace(x, d_map))
            d[k] = copy(l)
        else:
            d[k] = replace(d[k], d_map)
    return d


def test_interpolate():
    d = {'foo': 'hey',
         'bar': 'ho',
         'one': 'hey %(foo) ho %(bar)',
         'two': {'one': 'hey %(foo) ho %(bar)',
                 'two': 'raboof'}}
    d_interpolated = {'foo': 'hey',
                      'bar': 'ho',
                      'one': 'hey hey ho ho',
                      'two': {'one': 'hey hey ho ho',
                              'two': 'raboof'}}
    assert interpolate(d, d) == d_interpolated


def test_interpolate_int():
    d = {'foo': 1,
         'bar': 2,
         'one': 'hey %(foo) ho %(bar)',
         'two': {'one': 'hey %(foo) ho %(bar)',
                 'two': 'raboof'}}
    d_interpolated = {'foo': 1,
                      'bar': 2,
                      'one': 'hey 1 ho 2',
                      'two': {'one': 'hey 1 ho 2',
                              'two': 'raboof'}}
    assert interpolate(d, d) == d_interpolated


def test_interpolate_nested():
    d2 = {'modules': [{'fc': [{'source': '%(url_root)fc_optional.cmake'}]}], 'url_root': 'downloaded/downloaded_'}
    d2_interpolated = {'modules': [{'fc': [{'source': 'downloaded/downloaded_fc_optional.cmake'}]}], 'url_root': 'downloaded/downloaded_'}
    assert interpolate(d2, d2) == d2_interpolated
