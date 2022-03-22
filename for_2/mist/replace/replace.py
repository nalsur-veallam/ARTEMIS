from re import sub
from os.path import abspath, realpath, join, dirname

file = abspath(join(dirname(__file__), 'PARENT_MIST.txt'))
file_open = open(file, 'r')
file_read = file_open.read()
file_open.close()

new_file = abspath(join(dirname(__file__), 'MIST.txt'))
new_file_open = open(new_file, 'w')


def replace_content(dict_replace, target):
    """Based on dict, replaces key with the value on the target."""

    for check, replacer in list(dict_replace.items()):
        target = sub(check, replacer, target)
        # target = target.replace(check, replacer)

    return target


# check : replacer
dict_replace = {
    'd d': '33',
    'a d': '23',
    'a a': '22',
    'b a': '12',
    'b d': '13',
    'd a': '32',
    'a b': '21',
    'd b': '31',
    'b b': '11',
    'left to process: ': '',
    'MI: ': ''
}

new_content = replace_content(dict_replace, file_read)
new_file_open.write(new_content)
new_file_open.close()
