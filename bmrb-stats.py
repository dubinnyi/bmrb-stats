#!/usr/bin/python3 -u

#  Statistics from selected BMRB entries
#  Maxim A. Dubinnyi
#  11 March 2021

import pynmrstar
import argparse
import re
import sys
import os

LOG = sys.stdout
DBG = open(os.devnull, 'w')

def log_print(*args, **kwargs):
    kwargs['file'] = LOG
    print(*args, **kwargs)

def dbg_print(*args, **kwargs):
    print(*args, **kwargs, file=DBG)


def get_all_bmrb_entries(file = "/home/maxim/BMRB/all_entries.txt"):
    bmrb_id_list = []
    with open(file) as f:
        for line in f.readlines():
            bmrb_id_list.append(int(line.strip()))
    return bmrb_id_list

#
# regex for comama- and dash- separated selection of BMRB IDs, e.g.
# 30
# 30-40
# 30-40,50-60
# with spaces allowed in any place
text_to_expand_re = re.compile('(\d+(\s*-\s*\d+)?)(\s*,\s*\d+(\s*-\s*\d+)?)*')

def expand_id_text_to_list(id_text):
    id_list = []
    for text in id_text:
        if not text_to_expand_re.match(text):
            log_print("Wrong bmrb range specification: {}".format(text))
            continue
        else:
            for comma_separated in text.split(','):
                comma_separated = comma_separated.strip()
                dash_seperated = comma_separated.split('-')
                if len(dash_seperated) == 1:
                    id_list.append(int(dash_seperated[0]))
                elif len(dash_seperated) == 2:
                    id_first = int(dash_seperated[0])
                    id_last = int(dash_seperated[1])
                    id_range = range(id_first, id_last + 1)
                    log_print("Adding range {} - {}".format(id_first, id_last))
                    id_list.extend(id_range)
                else:
                    log_print("Wrong bmrb range specification: \'{}\' (\'{}\')".format(text, dash_seperated))
    return id_list


def select_bmrb_entries(all_entries_list, id_text_list):
    expanded_id_list = expand_id_text_to_list(id_text_list)
    all_entries_set = set(all_entries_list)
    expanded_id_set = set(expanded_id_list)
    selected_set = set.intersection(all_entries_set, expanded_id_set)
    return sorted(list(selected_set))


def get_entry_from_bmrb_id(bmrb_id):
    entry_file_name = "/home/maxim/BMRB/str/bmr{}_3.str".format(bmrb_id)
    return pynmrstar.Entry.from_file(entry_file_name)

polypeptide_re = re.compile('[Pp]olypeptide(L)')

def check_entity_is_polypeptide(entry, entity_id):
    entity_saveframes = entry.get_saveframes_by_category('entity')
    for saveframe in entity_saveframes:
        try:
            id = saveframe.tag_dict['id']
            type = saveframe.tag_dict['type']
            polytype = saveframe.tag_dict['polymer_type']
            if id == entity_id and type == 'polymer' and polypeptide_re.match(polytype):
                return True
        except:
            continue
    return False


lab_nucleus_uniform_re = re.compile('(U-)?(\d+%)?\s*(?P<nuc>(13C|15N|2H|17O))')

def parse_isotopic_labeling(lab_string, **kwargs):
    debug = True if 'debug' in kwargs else False
    func_name = sys._getframe().f_code.co_name
    if debug:
        print("debug {}, lab_string = {}".format(func_name, lab_string))
    lab_set = set()
    if lab_string[0] == '[':
        lab_string = lab_string[1:]
    if lab_string[-1] == ']':
        lab_string = lab_string[:-1]
    lab_splitted_sem = re.split('[,;]\s*',lab_string)
    if debug:
        print("debug {}, lab_splitted_sem = {}".format(func_name, lab_splitted_sem))
    for nuc_block in lab_splitted_sem:
        nuc_block = nuc_block.strip()
        add_nuc = '-'
        nuc_match = lab_nucleus_uniform_re.match(nuc_block)
        if nuc_match:
            add_nuc = nuc_match.group('nuc')
            lab_set.add(add_nuc)
        if debug:
            print("debug {},  {:>10} : {}".format(func_name, nuc_block, add_nuc))
    if debug:
        print("debug {}, lab_set = {}".format(func_name, lab_set))
    return lab_set


def print_entity_labeling(entry):
    entry_id = entry.get_tag('_Entry.ID')[0]
    entity_names = entry.get_tag('_Sample_component.Mol_common_name')
    entity_ids = entry.get_tag('_Sample_component.Entity_ID')
    entity_labeling = entry.get_tag('_Sample_component.Isotopic_labeling')
    entity_conc = entry.get_tag('_Sample_component.Concentration_val')
    entity_conc_units = entry.get_tag('_Sample_component.Concentration_val_units')
    flag_polypeptide = False
    lab_set = set()
    for name, id, lab, c, cu in \
            zip(entity_names, entity_ids, entity_labeling, entity_conc, entity_conc_units):
        if id.strip() != '.':
            if check_entity_is_polypeptide(entry, id):
                flag_polypeptide = True
                log_print("bmr{} LABEL: {} : {}, {} {}".format(entry_id, name, lab, c, cu))
                leb_entity_set =  parse_isotopic_labeling(lab)
                lab_set = lab_set.union(leb_entity_set)
    return flag_polypeptide, lab_set


shifts_data_re = re.compile('(?P<nuc>\d+[A-Z,a-z]+)\s+chemical shifts')

def print_shift_assignment(entry):
    shifts_nuc_set = set()
    entry_id = entry.get_tag('_Entry.ID')[0]
    datum_type = entry.get_tag('_Datum.Type')
    datum_count = entry.get_tag('_Datum.Count')
    for type,count in zip(datum_type, datum_count):
        log_print("bmr{} DATA : {} : {}".format(entry_id, type, count))
        shifts_nuc_match = shifts_data_re.match(type)
        if shifts_nuc_match:
            shifts_nuc_set.add(shifts_nuc_match.group('nuc'))
    return shifts_nuc_set


def assignment_strategy_heurystics(lab_nuc_set, shifts_nuc_set):

    labeled_assignment = []
    natural_assignment = []
    for nuc in ['1H', '15N', '13C']:
        if nuc in shifts_nuc_set:
            if nuc in lab_nuc_set:
                labeled_assignment.append(nuc)
            else:
                natural_assignment.append(nuc)
    ret_string = ''
    ret_string += "Labeled_"+"-".join(labeled_assignment) if labeled_assignment else ''
    if natural_assignment:
        ret_string += '_' if ret_string else ''
        ret_string += "Natural_" + "-".join(natural_assignment) if natural_assignment else ''
    if not ret_string:
        ret_string = 'NoAssignment'
    return ret_string


def print_labeling_and_assignment(entry):
    entry_id = entry.get_tag('_Entry.ID')[0]
    flag_polypeptide, lab_set = print_entity_labeling(entry)
    if flag_polypeptide:
        shifts_nuc_set = print_shift_assignment(entry)
        assignment_strategy = assignment_strategy_heurystics(lab_set, shifts_nuc_set)
        log_print("bmr{} TYPE : {}".format(entry_id, assignment_strategy))
    else:
        log_print("bmr{}: SKIP, not a polypeptide".format(entry_id))
        assignment_strategy = 'NotAPolypeptide'
    return assignment_strategy



def main():
    arg_parser = argparse.ArgumentParser(
        usage='Print labeling and assignment from selected BMRB entries')
    arg_parser.add_argument('idlist', type=str, help='bmrb IDs like 100-1000,12345', nargs='*')
    arg_parser.add_argument("--test", help="run simple test", action="store_true")
    arg_parser.add_argument("--log", help="print logging information for each record", action="store_true")
    args = arg_parser.parse_args()

    if not args.log:
        global LOG
        LOG=open(os.devnull,'w')
        print("suppress log_print")

    if args.test:
        test()
    else:
        all_entries_id = get_all_bmrb_entries()
        print("{} entries are in BMRB".format(len(all_entries_id)))

        if not args.idlist:
            selected_bmrb_ids = all_entries_id
        else:
            selected_bmrb_ids = select_bmrb_entries(all_entries_id, args.idlist)

        print("{} BMRB entries selected".format(len(selected_bmrb_ids)))
        strategy_counts = dict()
        strategy_ids = dict()
        for bmrb_id in selected_bmrb_ids:
            bmrb_entry = get_entry_from_bmrb_id(bmrb_id)
            strategy = print_labeling_and_assignment(bmrb_entry)
            if strategy in strategy_counts.keys():
                strategy_counts[strategy] += 1
                strategy_ids[strategy].append(bmrb_id)
            else:
                strategy_counts[strategy] = 1
                strategy_ids[strategy] = [bmrb_id]
            log_print('----')

        print("BMRB scan finished")
        for strategy in sorted(strategy_counts.keys()):
            print("{:>30} : {}".format(strategy, strategy_counts[strategy]))
            out_file = strategy + ".txt"
            with open(out_file, 'w') as keyf:
                for bmrb_id in strategy_ids[strategy]:
                    print(bmrb_id, file=keyf)


def test():
    test_lab_list = ["[U-100% 13C; U-100% 15N; 80% 2H]",
                     "[U-100% 15N; 80% 2H]",
                     "[80% 2H]"]
    for lab_u in test_lab_list:
        lab_set = parse_isotopic_labeling(lab_u, debug=True)
        print('---')

if __name__ == "__main__":
    main()
else:
    test()