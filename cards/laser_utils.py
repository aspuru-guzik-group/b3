# import os, sys, django
# os.environ["DJANGO_SETTINGS_MODULE"] = "djangochem.settings.postgres_docker"
# django.setup()

from pgmols.models import Mol, MolManager, Group
from rdkit.Chem import Draw
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw.rdMolDraw2D import MolDraw2DSVG
from rdkit.Chem.rdChemReactions import ReactionFromSmarts
from numpy.random import permutation
import xml.etree.ElementTree as ET
from django.conf import settings

from molgen.utils import put_smiles
import os
from datetime import datetime

INCHI_OPTION_STRING = " -RecMet  -FixedH "
SVG_DIR = getattr(settings, 'SVG_DIR')


# filter data from database
def filter_mols(group_name, *args, **kwargs):
    return Mol.objects.filter(group__name__exact=group_name, *args, **kwargs)

def delete_mols(mol_list, indexes = None):
    if indexes is None:
        mol_list.delete()
    else:
        for i in indexes:
            mol_list[i].delete()
        mol_list.update()

def list_projects():
    groups = Group.objects.all()
    return [g.name for g in groups]


def smiles_to_svg(smiles, **kwargs):
    img_size = kwargs.get('subImgSize', (300, 300))
    font_size = kwargs.get('font_size', 2.0)
    kekulize = kwargs.get('kekulize', True)

    mol = Chem.MolFromSmiles(smiles)
    mc = Chem.Mol(mol.ToBinary())
    if kekulize:
        try:
            Chem.Kekulize(mc)
        except:
            mc = Chem.Mol(mol.ToBinary())
    if not mc.GetNumConformers():
        rdDepictor.Compute2DCoords(mc)

    drawer = MolDraw2DSVG(img_size[0], img_size[1])
    drawer.SetFontSize(font_size)
    drawer.DrawMolecule(mc)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()

    # print(svg)
    ET.register_namespace("","http://www.w3.org/2000/svg")
    imgxml = ET.fromstring(svg)
    w = imgxml.attrib.pop('width')
    h = imgxml.attrib.pop('height')
    imgxml.set('viewBox', '0 0 {} {}'.format(w[:-2], h[:-2]))
    imgxml.set('xmlns', 'http://www.w3.org/2000/svg')
    imgxml.set('xmlns:xlink', 'http://www.w3.org/1999/xlink')
    svgstr = ET.tostring(imgxml).decode()
    return svgstr






def valid_smiles(smiles, project):
    # dupes
    inchikey = smiles_to_inchikey(smiles)
    if inchikey is None:
        return False

    mol = get_mols_inchikey(inchikey, project)[0]
    if mol is None:
        return True
    else:
        return False

def smiles_to_inchikey(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    tmp_inchikey = Chem.MolToInchi(mol, options=str(INCHI_OPTION_STRING))
    inchikey = Chem.InchiToInchiKey(tmp_inchikey)
    return inchikey

def get_mols_inchikey(inchikey, project):
    mols = Mol.objects.filter(group__name__exact=project, inchikey=inchikey)
    if mols.count() > 0:
        return mols
    else:
        return [None]
        

def full_path(key):
    return os.path.join(SVG_DIR, key+'.svg')


def fetch_svg(mol, **kwargs):
    svg_name = mol.inchikey + '.svg'
    svg_path = os.path.join(SVG_DIR, svg_name)

    if not os.path.isfile(svg_path):
        img = smiles_to_svg(mol.smiles, **kwargs)
        with open(svg_path, 'w') as f:
            f.write(img)
    else:
        with open(svg_path, 'r') as f:
            img = f.read()

    # return svg_name
    return img

    
def add_mol_to_db(smiles, project='tmp', info='', prop=''):
    return put_smiles(
        smiles, 
        project=project, 
        details={'info': info, "prop": prop}, 
        return_mol=True)


def has_all_elements(mol, match_elements, not_elements=[]):
    smi = mol.details['block'].replace('[Xe]', '').replace('[Be]', '').replace('[Ba]', '').lower()
    matches = [(el in smi) for el in match_elements] \
        + [(el not in smi) for el in not_elements]
    return sum(matches) == len(matches)



def run_reaction_svg(rxn_eqn, rxn_smarts, colors=[(.8,.9,1.),(1.,.9,.9)]):
    b_smi, x_smi = rxn_eqn.split('.')
    rxn = ReactionFromSmarts(rxn_smarts)

    # modify smiles
    def modify_smiles(smiles, remove, mapidx):
        smiles = smiles.replace(remove, '[*]')
        mol = Chem.MolFromSmiles(smiles)
        for a in mol.GetAtoms():
            a.SetProp('molAtomMapNumber', str(mapidx))
        return Chem.MolToSmiles(mol)
    b_smi = modify_smiles(b_smi, remove='[Xe]', mapidx=1)
    x_smi = modify_smiles(x_smi, remove='[Be]', mapidx=2)

    b_cnts = b_smi.count('Be')
    x_cnts = x_smi.count('Xe')


    b_mol = Chem.MolFromSmiles(b_smi)
    x_mol = Chem.MolFromSmiles(x_smi)

    # Error detection
    if b_cnts <= 0 or x_cnts <= 0:
        print(rxn_eqn)
        return 'ERROR'
    elif x_cnts == 1:
        p_mol = b_mol
        for _ in range(b_cnts):
            p_mol = rxn.RunReactants((p_mol, x_mol))[0][0]
    elif b_cnts == 1:
        p_mol = x_mol
        for _ in range(x_cnts):
            p_mol = rxn.RunReactants((b_mol, p_mol))[0][0]
    else:
        print('too much')
        return 'ERROR'

    # get p_smiles (repair mapping)
    for a in p_mol.GetAtoms():
        props = a.GetPropsAsDict()
        if props.get('old_mapno'):
            a.SetProp('molAtomMapNumber', a.GetProp('old_mapno'))
    p_smi = Chem.MolToSmiles(p_mol)

    # get svg
    svg = Draw.ReactionToImage(
        ReactionFromSmarts('{}.{}>>{}'.format(b_smi, x_smi, p_smi), useSmiles=True),
        highlightByReactant = True,
        highlightColorsReactants = colors,
        useSVG=True,)

    imgxml = ET.fromstring(svg.replace('svg:',''))
    w = imgxml.attrib.pop('width')
    h = imgxml.attrib.pop('height')
    imgxml.set('viewBox', '0 0 {} {}'.format(w[:-2], h[:-2]))
    imgxml.set('xmlns', 'http://www.w3.org/2000/svg')
    imgxml.set('xmlns:xlink', 'http://www.w3.org/1999/xlink')
    return ET.tostring(imgxml).decode()



def _fix_1_MIDA(smi):
    mol = Chem.MolFromSmiles(smi)
    midapatt = Chem.MolFromSmarts("B1(OC(=O)C[NR1](CC(=O)O1)C)")
    if mol.HasSubstructMatch(midapatt):
        while mol.HasSubstructMatch(midapatt):
            smarts="[B:1]~[O:2]~[C:3]~[C:4]~[N:5]>>[B-:1]1[*:2][*:3][*:4][N+:5]1"
            rxn = Chem.ReactionFromSmarts(smarts)
            mol = Chem.MolFromSmiles(smi)
            products = rxn.RunReactants((mol,))
            mol = products[0][0]
            smi = Chem.MolToSmiles(mol, isomericSmiles=True)
            mol = Chem.MolFromSmiles(smi)
    return smi

def fix_MIDA(smiles):
    smiles_list = smiles.split('.')
    for i, smi in enumerate(smiles_list):
        smiles_list[i] = _fix_1_MIDA(smi)
    return '.'.join(smiles_list)

def delete_svg(mol):
    svg_path = full_path(mol.inchikey)
    if os.path.isfile(svg_path):
        os.remove(svg_path)


### for adding smiles! add mol through smiles
from . import laser_search as search
def getmol(smiles): return Chem.MolFromSmiles(smiles)
def getsmiles(mol): return Chem.MolToSmiles(mol)
trusted_vendors = [
    'Sigma-Aldrich', 'Combi-Blocks', 'Oakwood Products', 
    'Enamine', 'Alfa Chemistry', 'Luminescence Technology Corp. (Lumtec)']
b_frags = [
    '[H]B1Oc2ccccc2O1',
    '[H]B1OC(=O)CN(C)CC(=O)O1',
    '[H][B-]12OC(=O)C[N+]1(C)CC(=O)O2',
    '[H]B1OC(C)(C)C(C)(C)O1',
    '[H]B1Nc2cccc3cccc(c23)N1',
    '[H]B(O)O',
]
b_patts = [
    "B1(OC(=O)C[NR1](CC(=O)O1)C)", 
    "[B-]12(OC(=O)C[N+]2(CC(=O)O1)C)",
    "B(O)O",
]
def replace_pattern(smiles, patt, repl, replaceAll = False):
    mol = getmol(smiles)
    rmol = Chem.ReplaceSubstructs(mol, getmol(patt), getmol(repl), replaceAll=replaceAll)
    try:
        smiles_list = [getsmiles(Chem.RemoveHs(rm)) for rm in rmol]
    except:
        smiles_list = [smiles]
    return smiles_list
def x_end(s):
    s = s.replace('Br', '[Xe]').replace('I', '[Xe]')
    return s
def b_end(smi):

    mol = Chem.MolFromSmiles(smi)
    for b_patt in b_patts:
        b_mol = Chem.MolFromSmarts(b_patt)
        while mol.HasSubstructMatch(b_mol):
            smarts="[B:1][C,c:2]>>[Be]/[*:2].[*:1]"
            rxn = Chem.ReactionFromSmarts(smarts)
            mol = Chem.MolFromSmiles(smi)
            products = rxn.RunReactants((mol,))
            mol = products[0][0]
            smi = Chem.MolToSmiles(mol, isomericSmiles=True)
            smi = smi.replace("[BeH]","[Be]")
            mol = Chem.MolFromSmiles(smi)
    return smi


def get_info_by_smiles(smiles):
    cid = search.get_cid_by_smiles(smiles)
    if cid:
        cid = str(cid)
    else:
        return None, None, None, [None]
    block = b_end(x_end(smiles))
    tag = ''
    for flag in ['Be', 'Ba', 'Xe']:
        if flag in block:
            tag += flag
    if tag == '':
        tag = 'None'
    urls = search.get_vendors_from_cid(cid, trusted_vendors)
    return smiles, cid, urls, [tag]




def min_price(a, b):
    ap = 1e20 if a['price'] == 'Nan' else a['price']
    bp = 1e20 if b['price'] == 'Nan' else b['price']
    if ap < bp:
        return a
    else:
        return b
def update_prices(mol, new_prices):
    for i, price in enumerate(mol.details['prices']):
        mol.details['prices'][i] = min_price(price, new_prices[i])

def add_block_by_smiles(smiles, project):
    "returns inchikey, message"
    try:
        smiles = getsmiles(getmol(smiles))
    except:
        print('Invalid smiles')
        return None, 'Invalid smiles'
    
    # check if molecules already exist
    db_mols = filter_mols(project)
    block = getsmiles(getmol(b_end(x_end(smiles))))
    match_mol = None
    for mol in db_mols:
        if mol.details['block'] == block:
            match_mol = mol
            break
            
    _, cid, urls, tags = get_info_by_smiles(smiles)
    smiles = fix_MIDA(smiles)
    if urls in ['', None]:
        print('No trusted vendors')
        return None, 'No trusted vendors'
    if match_mol is not None:
        if cid not in match_mol.details['cid']:
            prices = search.get_vendor_prices(urls)
            match_mol.details['smiles'] += ('.'+smiles)
            match_mol.smiles += ('.'+smiles)
            match_mol.details['cid'] += ('!!'+cid)
            match_mol.details['url'] += ('!!'+urls)
            update_prices(match_mol, prices)
            match_mol.save()
            svg_path = full_path(match_mol.inchikey)
            if os.path.isfile(svg_path):
                os.remove(svg_path)
        return match_mol.inchikey, 'Exist'
            
    else:
        prices = search.get_vendor_prices(urls)
        details={
            'block': block,
            'smiles': smiles,
            'cid': cid, 
            'url': urls, 
            'prices': prices,
            "prop": 'U',
        }
        new_mol = put_smiles(
            block+'.'+smiles, 
            project = project,
            tags = tags,
            details = details,
            return_mol=True,
        )
        return new_mol['mol'].inchikey, 'Success'



def add_block_by_cas(cas, project):
    "returns list of inchikey, message"
    cas = cas.strip().replace(' ', '')
    smiles = search.get_smiles_from_cas(cas)

    if smiles:
        return add_block_by_smiles(smiles, project)
    else:
        return None, 'Invalid CAS'




# for views
def now():
    return f'({datetime.now():%Y/%m/%d %H:%M:%S})'

def get_ip(request):
    x_forwarded_for = request.META.get('HTTP_X_FORWARDED_FOR')
    if x_forwarded_for:
        ip = x_forwarded_for.split(',')[0]
    else:
        ip = request.META.get('REMOTE_ADDR')
    return ip

def printlog(viewfunc, request):
    print(f'{now()} [{get_ip(request)}] {viewfunc} {request.method}')

