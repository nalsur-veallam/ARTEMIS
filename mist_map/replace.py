from re import sub
from os.path import abspath, realpath, join, dirname

file = abspath(join(dirname(__file__), 'PARENT_MIE.txt'))
file_open = open(file, 'r')
file_read = file_open.read()
file_open.close()

new_file = abspath(join(dirname(__file__), 'MIE.txt'))
new_file_open = open(new_file, 'w')


def replace_content(dict_replace, target):
    """Based on dict, replaces key with the value on the target."""

    for check, replacer in list(dict_replace.items()):
        target = sub(check, replacer, target)
        # target = target.replace(check, replacer)

    return target


# check : replacer
dict_replace = {
    'dihedral-dihedral': '33',
    'angle-dihedral': '23',
    'angle-angle': '22',
    'bond-angle': '12',
    'bond-dihedral': '13',
    'bond-bond': '11'
}

new_content = replace_content(dict_replace, file_read)
new_file_open.write(new_content)
new_file_open.close()
