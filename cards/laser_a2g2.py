import django
django.setup()
from pgmols.models import Mol, Group
from molgen.utils import put_mol, put_smiles
# from .laser_utils import getmol, getsmiles



# filter data from database
# tags__contains
def filter_mols(project=None, *args, **kwargs):
    if project is not None:
        return Mol.objects.filter(group__name__exact=project, *args, **kwargs)
    else:
        return Mol.objects.all()

def delete_mols(mol_list, indexes = None):
    if indexes is None:
        mol_list.delete()
    else:
        for i in indexes:
            mol_list[i].delete()
        mol_list.update()

def delete_project(project):
    "this only deltes all the molecules in the project"
    mol_list = filter_mols(project)
    delete_mols(mol_list)

def list_projects():
    groups = Group.objects.all()
    return [g.name for g in groups]

def copy_project(old_project, new_project):
    old_mols = filter_mols(old_project)
    for mol in old_mols:
        details = mol.details.copy()
        details.pop('mass')
        details.pop('stoichiometry')
        details.pop('electron_count')
        details.pop('molecular_charge')

        tags = mol.tags.copy()
        put_smiles(mol.smiles, project=new_project, tags=tags, details=details)
    return

def add_mols(mol_list, project, **kwargs):
    for mol in mol_list:
        put_mol(mol, project = project, **kwargs)
    return

def modify_mol_details(mol, key, func, update=False):
    if mol.details.get(key, None) is None or update:
        mol.details[key] = func(mol)
        mol.save()
    else:
        print('Did not update mol.details["{}"]'.format(key))

def add_smiles_list(smiles_list, project, **kwargs):
    for s in smiles_list:
        put_smiles(s, project = project, **kwargs)
    return

