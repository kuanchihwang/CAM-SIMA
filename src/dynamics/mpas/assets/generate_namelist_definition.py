#!/usr/bin/env python3

'''
Generate XML namelist definition file for MPAS dynamical core in CAM-SIMA.
'''

import argparse
import textwrap
import xml.etree.ElementTree as ET

EXCLUDED_NAMELIST_GROUP = [
    'limited_area',
    'physics'
]
EXCLUDED_NAMELIST_OPTION = [
    'config_calendar_type',
    'config_do_restart',
    'config_run_duration',
    'config_start_time',
    'config_stop_time'
]
INDENT_PER_LEVEL = ' ' * 4
NEW_PREFIX = 'mpas_'
OLD_PREFIX = 'config_'

def parse_argument() -> argparse.Namespace:
    '''
    Parse command line arguments.
    '''

    parser = argparse.ArgumentParser(
        description='Generate XML namelist definition file for MPAS dynamical core in CAM-SIMA.'
    )

    parser.add_argument(
        '-r', '--registry',
        default='Registry.xml',
        type=str,
        required=False,
        help='XML MPAS registry file.',
        dest='reg_xml'
    )
    parser.add_argument(
        '-n', '--namelist',
        default='Namelist.xml',
        type=str,
        required=False,
        help='XML namelist definition file.',
        dest='nml_xml'
    )
    parser.add_argument(
        '-s', '--schema',
        default=None,
        type=str,
        required=False,
        help='XML schema for namelist definition file.',
        dest='nml_xsd'
    )

    argument = parser.parse_args()

    return argument

def parse_xml(xml_file: str) -> ET.ElementTree:
    '''
    Parse XML file into element tree.
    '''

    xml_et = ET.parse(xml_file)

    return xml_et

def validate_xml(xml_file: str, xsd_file: str) -> bool:
    '''
    Validate XML file against XSD file.
    '''

    # Only import `xmlschema` if XML validation is requested. Silence pylint about it.
    import xmlschema # pylint: disable=import-outside-toplevel

    xml_schema = xmlschema.XMLSchema(xsd_file)

    return xml_schema.is_valid(xml_file)

def transform_name(name: str) -> str:
    '''
    Change prefix of namelist option/group name.
    '''

    while name.startswith(OLD_PREFIX):
        name = name[len(OLD_PREFIX):]

    while name.startswith(NEW_PREFIX):
        name = name[len(NEW_PREFIX):]

    name = NEW_PREFIX + name

    return name

def translate_element_tree(reg_xml_et: ET.ElementTree) -> ET.ElementTree:
    '''
    Translate MPAS registry into namelist definition.
    '''

    # `entry_id_pg` is the root element in namelist definition.
    entry_id_pg = ET.Element('entry_id_pg', {'version': '0.1'})

    comment = ET.Comment(
        '\n' +
        INDENT_PER_LEVEL * 2 + 'MPAS dycore' + '\n' +
        '\n' +
        INDENT_PER_LEVEL * 2 + 'Notes to developers/maintainers:' + '\n' +
        INDENT_PER_LEVEL * 2 + 'This file is auto-generated from MPAS registry. Do not edit directly.' + '\n' +
        INDENT_PER_LEVEL * 2 + 'Instead, use the Python script at `src/dynamics/mpas/assets/generate_namelist_definition.py`.' + '\n' +
        INDENT_PER_LEVEL
    )
    entry_id_pg.append(comment)

    for namelist_group in reg_xml_et.findall('nml_record'):
        if namelist_group.attrib['name'].strip().lower() in EXCLUDED_NAMELIST_GROUP:
            continue

        for namelist_option in namelist_group.findall('nml_option'):
            if namelist_option.attrib['name'].strip().lower() in EXCLUDED_NAMELIST_OPTION:
                continue

            # The `entry_id_pg` root element contains many `entry` elements.
            # Each `entry` element describes a namelist option, indicated by its `id` attribute.
            entry = ET.SubElement(entry_id_pg, 'entry', {'id': transform_name(namelist_option.attrib['name'].strip().lower())})

            # The `category` element.
            sub_element = ET.SubElement(entry, 'category')
            sub_element.text = 'mpas'

            # The `desc` element.
            desc_text = ' '.join(namelist_option.attrib['description'].strip('.').split())
            desc_text = desc_text[0].upper() + desc_text[1:]
            desc_text = '\n' + textwrap.fill(desc_text, 80, initial_indent=INDENT_PER_LEVEL * 3, subsequent_indent=INDENT_PER_LEVEL * 3) + '\n' + INDENT_PER_LEVEL * 2

            sub_element = ET.SubElement(entry, 'desc')
            sub_element.text = desc_text

            # The `group` element.
            sub_element = ET.SubElement(entry, 'group')
            sub_element.text = transform_name(namelist_group.attrib['name'].strip().lower())

            # The `type` element.
            # The `values` element and its containing `value` element.
            type_text = namelist_option.attrib['type'].strip().lower()

            if type_text.startswith('c'):
                type_text = 'char*256'
                value_text = namelist_option.attrib['default_value'].strip()
            elif type_text.startswith('i'):
                type_text = 'integer'
                value_text = namelist_option.attrib['default_value'].strip().lower()
            elif type_text.startswith('l'):
                type_text = 'logical'
                value_text = namelist_option.attrib['default_value'].strip().lower()

                if value_text.startswith(('t', '.t')):
                    value_text = '.true.'
                elif value_text.startswith(('f', '.f')):
                    value_text = '.false.'
                else:
                    raise ValueError('Invalid value')
            elif type_text.startswith('r'):
                type_text = 'real'
                value_text = namelist_option.attrib['default_value'].strip().lower()

                if value_text.startswith('.'):
                    value_text = '0' + value_text

                if value_text.endswith('.'):
                    value_text = value_text + '0'

                i = value_text.find('.d')

                if i != -1:
                    value_text = value_text[:i + 1] + '0' + value_text[i + 1:]

                i = value_text.find('.e')

                if i != -1:
                    value_text = value_text[:i + 1] + '0' + value_text[i + 1:]
            else:
                raise ValueError('Invalid type')

            sub_element = ET.SubElement(entry, 'type')
            sub_element.text = type_text

            sub_element = ET.SubElement(entry, 'values')
            sub_element = ET.SubElement(sub_element, 'value')
            sub_element.text = value_text

    # Sort the `entry` elements for result stability except for the comment element at index 0.
    entry_id_pg[1:] = sorted(entry_id_pg.iterfind('entry'), key=lambda entry: entry.get('id'))
    nml_xml_et = ET.ElementTree(entry_id_pg)

    return nml_xml_et

def write_element_tree(xml_file: str, xml_et: ET.ElementTree) -> None:
    '''
    Write element tree into XML file.
    '''

    ET.indent(xml_et, space=INDENT_PER_LEVEL)

    xml_et.write(
        xml_file,
        encoding='UTF-8',
        xml_declaration=True,
        default_namespace=None,
        method='xml',
        short_empty_elements=False
    )

    # The file written by `ElementTree.write()` contains no newline at end of file. Add it manually.
    with open(xml_file, 'a', encoding='utf-8') as nml_xml_file:
        nml_xml_file.write('\n')

if __name__ == '__main__':
    arg = parse_argument()
    reg_xml_element_tree = parse_xml(arg.reg_xml)
    nml_xml_element_tree = translate_element_tree(reg_xml_element_tree)

    write_element_tree(arg.nml_xml, nml_xml_element_tree)
    print('Generated ' + arg.nml_xml)

    if arg.nml_xsd is not None:
        if validate_xml(arg.nml_xml, arg.nml_xsd):
            print('Successfully validated ' + arg.nml_xml + ' against ' + arg.nml_xsd)
        else:
            print('Failed to validate ' + arg.nml_xml + ' against ' + arg.nml_xsd)
