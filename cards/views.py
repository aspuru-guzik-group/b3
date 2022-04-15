from django.http import HttpResponse, HttpResponseRedirect, JsonResponse
from django.contrib.auth import authenticate, login
from django.db.models import Q
# from .laser_utils import filter_mols, fetch_svg, valid_smiles, add_mol_to_db, smiles_to_inchikey
from . import laser_utils as ut
from django.shortcuts import render, redirect
from .forms import MolForm, BlockForm
from django.contrib.auth.decorators import login_required
import csv

from numpy.random import choice


# Usage: project
# Usage: project/<project_name>
# Usage: project/<project_name>/data.csv
# Usage: project/<project_name>/new

pubchem_url = 'https://pubchem.ncbi.nlm.nih.gov/compound/{}'
block_style = {
    'subImgSize': (400, 500),
    'font_size': 6
}

@login_required
def projects(request):
    projects = ut.list_projects()
    return render(request, 'projects.html', { 
        'projects': projects,
    })

@login_required
def projectpage(request, project=None):
    ut.printlog(project, request)

    db_mols = ut.filter_mols(project)
    if request.method == 'POST':
        inchikey = request.POST.get('inchikey')
        flag = request.POST.get('flag')

        mols = ut.get_mols_inchikey(inchikey, project=project)
        if mols[0] is not None:
            for m in mols:
                print('Selection update:', inchikey, flag)
                m.details.update(prop=flag)
                m.save()

    # extra fitlers
    tag = request.GET.get('tag', '')
    flag = request.GET.get('flag', None)
    elements = request.GET.get('elements', None)
    perrow = int(request.GET.get('perrow', 8))
    inchikeys = request.GET.get('inchikey', '')
    molrange = request.GET.get('range', None)
    print('Get params:', [tag, flag, perrow, inchikeys, elements, molrange])

    if tag != '':
        db_mols = [mol for mol in db_mols if tag in mol.tags]
    if flag is not None:
        flags = flag.split(',')
        db_mols = [mol for mol in db_mols if mol.details['prop'] in flags]
    if elements is not None:
        elements = list(map(str.lower, elements.split(',')))
        match_elements = [el for el in elements if el[0]!='^']
        not_elements = [el[1:] for el in elements if el[0]=='^']
        db_mols = [mol for mol in db_mols if ut.has_all_elements(mol, match_elements, not_elements)]
    if molrange is not None:
        ranges = molrange.split('-')
        if len(ranges) == 2:
            db_mols = db_mols[int(ranges[0]):int(ranges[1])]

    margin = 0.2
    cardwidth = round(100./perrow - 2*margin, 5)
    infowidth = 2*cardwidth+2*margin
    smiles = lambda mol: ''

    if inchikeys != '':
        db_mols = ut.filter_mols(project, inchikey__in=inchikeys.split(','))
        if len(db_mols) == 1:
            cardwidth = 25-2*margin
            smiles = lambda mol: mol.details.get('smiles', '') # list smiles in the card


    get_links = lambda cids: [
        {'name': v, 'url': pubchem_url.format(v)} for i, v in enumerate(cids)]

    # for url tags
    params = []
    if tag:
        params.append('tag=' + tag)
    if flag:
        params.append('flag=' + flag)
    if inchikeys:
        if len(db_mols) <= 5:
            params.append('inchikey=' + inchikeys)
    if len(params) > 0:
        download_params = '?' + '&'.join(params)
    else:
        download_params = ''

    # define molecule data to pass
    def mol_context(mol):
        return {
            'imgsrc': ut.fetch_svg(mol, **block_style),
            'inchikey': mol.inchikey,
            'pubchem': get_links(mol.details['cid'].split('!!')),
            'prices': mol.details.get('prices', []),
            'smiles': smiles(mol),
            'f': 'active' if mol.details.get('prop', '')=='F' else '',
            'y': 'active' if mol.details.get('prop', '')=='Y' else '',
            'n': 'active' if mol.details.get('prop', '')=='N' else '',
            'u': 'active' if mol.details.get('prop', '')=='U' else '',
        }

    # get molecules
    mols = map(mol_context, db_mols)
    return render(request, 'projectmols.html', { 
        'molecules': mols, 
        'cardwidth': cardwidth, 
        'infowidth': infowidth,
        'margin': margin,
        'params': download_params,
        'project': project,
    })

@login_required
def projectpage_download(request, project=None):
    ut.printlog(f'{project}_download', request)

    # filter data
    db_mols = ut.filter_mols(project)
    inchikeys = request.GET.get('inchikeys', '')
    inchikey = request.GET.get('inchikey', None)
    # single molecule fetch
    if inchikeys != '':
        inchikeys = inchikeys.split(',')
        db_mols = [mol for mol in db_mols if mol.inchikey in inchikeys]

    elif inchikey is not None:
        db_mols = [mol for mol in db_mols if mol.inchikey == inchikey]

    # other flag fetch
    else:
        tag = request.GET.get('tag', '')
        flag = request.GET.get('flag', None)
        if tag != '':
            db_mols = [mol for mol in db_mols if tag in mol.tags]
        if flag is not None:
            flags = flag.split(',')
            db_mols = [mol for mol in db_mols if mol.details['prop'] in flags]

    response = HttpResponse(content_type='text/csv')
    response['Content-Disposition'] = 'attachment; filename="data.csv"'
    writer = csv.writer(response)

    writer.writerow([
        'inchikey', 'block', 'smiles', 'flag', 
        'lowest_price', 'pubchem_cid', 'urls',
    ])

    for mol in db_mols:
        lowest_price = 1e9
        for prices in mol.details.get('prices', [{'price': 'Nan'}]):
            if prices['price'] != 'Nan' and prices['price'] < lowest_price:
                lowest_price = prices['price']
        lowest_price = lowest_price if lowest_price < 9e8 else ''

        writer.writerow([
            mol.inchikey, 
            mol.details['block'], 
            mol.details['smiles'], 
            mol.details['prop'], 
            lowest_price, 
            mol.details['cid'], 
            mol.details['url'], 
        ])

    return response


@login_required
def projectpage_addpage(request, project=None):
    message, list_message, cas_message = '', '', ''
    ut.printlog(f'{project}_addpage', request)

    if request.method == 'POST':
        form = BlockForm(request.POST)
        if form.is_valid():
            if form.cleaned_data['button_type'] == 'add':
                smiles = form.cleaned_data['smiles']
                print('adding:', smiles)
                inchikey, message = ut.add_block_by_smiles(smiles, project)
                if inchikey:
                    return HttpResponseRedirect(f'/cards/project/{project}?inchikey='+inchikey)

            elif form.cleaned_data['button_type'] == 'list': # add a smiles list
                inchikeys = []
                smiles_list = form.cleaned_data['smiles_list']
                for smiles in smiles_list.split('\n'):
                    smiles = smiles.strip()
                    if smiles == '':
                        continue
                    print('adding:', smiles)
                    inchikey, _ = ut.add_block_by_smiles(smiles, project)
                    if inchikey:
                        inchikeys.append(inchikey)
                if len(inchikeys) > 0:
                    return HttpResponseRedirect(f'/cards/project/{project}?inchikey='+','.join(inchikeys))
                else:
                    list_message = 'Empty list'

            elif form.cleaned_data['button_type'] == 'cas': # add a smiles list
                inchikeys = []
                cas_list = form.cleaned_data['cas_list']

                for cas in cas_list.split('\n'):
                    cas = cas.strip().replace(' ', '')
                    print('adding:', cas.strip())
                    inchikey, _ = ut.add_block_by_cas(cas, project)
                    if inchikey:
                        inchikeys.append(inchikey)

                if len(inchikeys) > 0:
                    return HttpResponseRedirect(f'/cards/project/{project}?inchikey='+','.join(inchikeys))
                else:
                    cas_message = 'Empty list'

            else:
                raise ValueError('wrong button values?')

        form = BlockForm()
        return render(request, 'blocknew.html', {
                'form': form, 
                'message': message, 
                'list_message': list_message, 
                'cas_message': cas_message
        })

    else:
        form = BlockForm()
        return render(request, 'blocknew.html', 
            {'form': form, 'message': '', 'list_message': '', 'cas_message': ''})





