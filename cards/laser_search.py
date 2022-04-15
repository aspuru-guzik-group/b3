import re, time
from datetime import datetime, date
import urllib.request
import urllib.parse
import json, time
import xml.etree.ElementTree as ET
from io import StringIO
import pandas as pd
import numpy as np
from selenium import webdriver
from selenium.webdriver.common.keys import Keys


url_substructure = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/substructure/smiles/{}/JSON'
url_smiles = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{}/JSON'
url_smiles_post = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/JSON'
url_fastsubstructure = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsubstructure/smiles/{}/property/CanonicalSMILES/JSON?StripHydrogen=true'
url_cids = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/listkey/{}/cids/JSON'
url_cids_post = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/property/CanonicalSMILES/TXT'

## substructure search ##
def url_request(url, post=None, fmt=None, headers = {}, timeout=600):
    req =  urllib.request.Request(url, headers=headers)
    post = post.encode() if post is not None else None
    # http request
    try:
        resp = urllib.request.urlopen(req, data=post, timeout=timeout).read().decode('utf-8')
    except urllib.error.HTTPError as e:
        resp = e.read()
        if fmt is not None and fmt.upper() == 'JSON':
            return json.loads(resp)
        else:
            return resp
    except UnicodeDecodeError:
        resp = urllib.request.urlopen(req, data=post, timeout=timeout).read().decode('ISO-8859-1')

    # format response
    if fmt is None:
        return resp
    elif fmt.upper() == 'JSON':
        return json.loads(resp)
    elif fmt.upper() == 'CSV':
        csv_io = StringIO(resp)
        return pd.read_csv(csv_io)
    elif fmt.upper() == 'TXT':
        return [s for s in resp.split('\n') if len(s)>0]
    elif fmt.upper() == 'XML':
        return ET.fromstring(resp)
    else:
        raise TypeError('fmt has wrong data type')

def get_mols_by_substructure(sub):
    sub = sub.replace('#', '%23')
    response = url_request(url_fastsubstructure.format(sub), fmt='json')
    if response.get('Fault', None):
        print(' '+response['Fault']['Message'], '...', end='')
        return pd.DataFrame(columns=['cid', 'smiles'])
    else:
        cid_smiles = pd.DataFrame(response['PropertyTable']['Properties'])
        cid_smiles.rename(columns={'CID': 'cid', 'CanonicalSMILES': 'smiles'}, inplace = True)
        return cid_smiles


def get_cid_by_smiles(smiles, message=False):
    response = url_request(url_smiles_post, post='smiles='+smiles, fmt='json')
    try:
        cid = response['PC_Compounds'][0]['id']['id']['cid']
        return cid
    except:
        return None




def get_cids_by_substructure(sub, message=False):
    response = url_request(url_substructure.format(sub).replace('#', '%23'), fmt='json')
    key = response['Waiting']['ListKey']
    time.sleep(0.5)

    waiting = True
    while waiting:
        response = url_request(url_cids.format(key), fmt='json')
        waiting = response.get('Waiting', False)
        if message:
            print(waiting)
        time.sleep(1)
    if response.get('Fault', False):
        if message:
            print(response['Fault']['Message'], '...', end='')
    id_list = response.get('IdentifierList', {'CID': []})['CID']
    return id_list

def get_smiles_from_cids(cids):
    cids_post = 'cid=' + ','.join([str(c) for c in cids])
    return url_request(url_cids_post, post=cids_post, fmt='txt')

## download smiles with cids ##
def fragments_by_substructures(pair_sub):
    # substructure search
    sub_dict = {}
    mols_df = pd.DataFrame(columns=['cid', 'smiles'])

    for e, e_sub in pair_sub:
        e_cid_smiles = pd.DataFrame(columns=['cid', 'smiles'])
        sub_dict[e] = []
        for sub in e_sub:
            print('matching', sub, '...', end='')
            cid_smiles = get_mols_by_substructure(sub)
            mols_df = mols_df.append(cid_smiles)
            print(' found', cid_smiles.shape[0], 'mols.')
            
            sub_dict[e].append((sub, cid_smiles['cid'].values))
    
    # remove duplicates and sort by cid
    mols_df.sort_values(by='cid', inplace=True)
    mols_df.drop_duplicates(subset='cid', inplace=True)
    mols_df.set_index('cid', inplace=True)
    
    mols_df = mols_df.assign(frags='')
    for e in sub_dict.keys():
        sub_cids_pairs = sub_dict[e]
        mols_df = mols_df.assign(**{e: 0})
        for sub, cids in sub_cids_pairs:
            mols_df.loc[cids, 'frags'] += (sub+'.')
            mols_df.loc[cids, e] = 1
    return mols_df


## download molecules from esearch ##
def save_mols_from_esearch(fname, term, info_func=None):
    batch_size = 100000
    vendor_url = 'https://eutils.ncbi.nlm.nih.gov/eutils/esearch.fcgi?db=pccompound&usehistory=y'

    fetch_url = 'https://eutils.ncbi.nlm.nih.gov/eutils/efetch.fcgi?db=pccompound&rettype=uilist&retmode=tet&query_key={key}&retstart={i}&retmax={batch_size}&WebEnv={webenv}'
    
    # get the search history
    et = url_request(vendor_url, post='term={}'.format(term), fmt='xml')
    webenv = list(et.iter('WebEnv'))[0].text
    count = int(list(et.iter('Count'))[0].text)
    qkey = list(et.iter('QueryKey'))[0].text
    
    # fetch the history one batch at a time due to its limits
    header = True
    nums = 0
    with open(fname, 'a') as f:
        while nums < count:
            batch_url = fetch_url.format(i=nums, webenv=webenv, batch_size=batch_size, key=qkey)
            resp = url_request(batch_url)

            try:
                cids = np.fromstring(resp, dtype='u4', sep='\n')
                nums += len(cids)
                smiles_list = get_smiles_from_cids(cids)
                df = pd.DataFrame({'cid': cids, 'smiles': smiles_list})
                if info_func is not None:
                    df_info = df.apply(info_func, axis=1)
                    df = pd.concat([df, df_info], axis=1) # .drop(columns=['Unnamed: 0'])
                df.to_csv(f, header=header, index=False)
                header=False

                print('\rFinished {}/{}, {:d}%'.format(nums, count, int(100.*nums/count)), end='')
                time.sleep(0.2)
            except:
                print('\nCannot read {} ~ {}, skipping\n'.format(nums, nums+batch_size))
                nums += batch_size


def get_cids_from_esearch(term):
    batch_size = 100000
    vendor_url = 'https://eutils.ncbi.nlm.nih.gov/eutils/esearch.fcgi?db=pccompound&usehistory=y&term={}'
    fetch_url = 'https://eutils.ncbi.nlm.nih.gov/eutils/efetch.fcgi?db=pccompound&rettype=uilist&retmode=text&query_key={key}&retstart={i}&retmax={batch_size}&WebEnv={webenv}'
    
    # get the search history
    et = url_request(vendor_url.format(term), fmt='xml')
    webenv = list(et.iter('WebEnv'))[0].text
    count = int(list(et.iter('Count'))[0].text)
    qkey = list(et.iter('QueryKey'))[0].text
    
    # fetch the history one batch at a time due to its limits
    cids_from_search = np.zeros((0,), dtype='u4')
    nums = 0
    while nums < count:
        batch_url = fetch_url.format(i=nums, webenv=webenv, batch_size=batch_size, key=qkey)
        resp = url_request(batch_url)
        try:
            cids = np.fromstring(resp, dtype='u4', sep='\n')
            nums += len(cids)
            cids_from_search = np.r_[cids_from_search, cids]
            print('\rFinished {}/{}, {:d}%'.format(nums, count, int(100.*nums/count)), end='')
            time.sleep(0.2)
        except:
            print('\nCannot read {} ~ {}, skipping\n'.format(nums, nums+batch_size))
            nums += batch_size
    return cids_from_search[::-1]

TRUSTED_VENDORS = ['Sigma-Aldrich', 'Combi-Blocks', 'Oakwood Products', 'Enamine', 'Alfa Chemistry', 'Luminescence Technology Corp. (Lumtec)']
def get_vendors_from_cid(cid, trusted_vendors=TRUSTED_VENDORS):
    vendor_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/categories/compound/{}/JSON'
    res = url_request(vendor_url.format(cid), fmt='json')
    source_urls = []

    if res.get("SourceCategories", None) is None:
        return '$$'.join(source_urls)

    for cat in res['SourceCategories']['Categories']:
        if cat['Category'] == 'Chemical Vendors':
            for source in cat['Sources']:
                if source['SourceName'] in trusted_vendors:
                    url = source.get('SourceRecordURL', None)
                    if url is not None:
                        source_urls.append(url)
    return '$$'.join(source_urls)

def get_smiles_from_cas(cas):
    # cas_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/cids/JSON'
    cas_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/property/CanonicalSMILES/JSON'
    data = url_request(cas_url.format(cas), fmt='JSON')

    if 'PropertyTable' in data:
        return data['PropertyTable']['Properties'][0]['CanonicalSMILES']
    else:
        return None


def has_trusted_vendors(cids):
    batch_size = 100000
    return

# version 2
cids_url_post = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/xrefs/SBURL/JSON'
vendor_re = re.compile(r'\w+\://[-a-zA-Z]+\.([-a-zA-Z]+)\.[-a-zA-Z]+.+')

selected_vendors = [
    'sigmaaldrich', 'enaminestore', 'alfa.com',
    'alfa-chemistry', 'fishersci', 'combi-blocks',
    'lumtec', 'oakwoodchemical', 'tcichemicals',
    'avantorinc', ]
def has_selected_vendors(url):
    for vendor_str in selected_vendors:
        if vendor_str in url:
            return True
    return False
        
def get_batched_vendors(cids):
    cids_post = 'cid=' + ','.join([str(c) for c in cids])
    data = url_request(cids_url_post, post=cids_post, fmt='JSON')
    urls = []
    for d in data['InformationList']['Information']:
        urls.append('$$'.join(filter(has_selected_vendors, d.get('SBURL', []))))
    return urls

## post matching and filtering
def filter_candidates(fname, funcs, dbname='data/frags_candidates.csv', chunksize=10000):
    with open(dbname, 'r') as f:
        data = f.readlines()
        mol_counts = len(data) - 1

    header = True
    counts = 0
    with open(fname, 'a') as f:
        for chunk in pd.read_csv(dbname, chunksize=chunksize):
            counts+= chunk.shape[0]
            chunk_info = chunk.apply(funcs, axis=1)
            columns = list(np.sort(chunk.columns.values)) + list(np.sort(chunk_info.columns.values))
            if chunk_info.shape[1] > 0:
                chunk = pd.concat([chunk.loc[chunk_info.index.values], chunk_info], axis=1).dropna()
                chunk[columns].to_csv(f, header=header, index=False)
                header = False
            print('\rUpdated {}/{} data, {:d}%     '.format(counts, mol_counts, int(100.*counts/mol_counts)), end='')


def filter_data(df, filters):
    if df.shape[0] <=0:
        return df
    for filt in filters:
        df = df[filt(df)]
    return df

def data_operation(df, operations):
    if df.shape[0] <= 0:
        return df
    for oper in operations:
        df = oper(df)
    return df


## get prices from vendors
sigma_aldrich_url = 'https://www.sigmaaldrich.com/catalog/PricingAvailability.do?productNumber={}&brandKey=ALDRICH'
oakwood_url = 'http://www.oakwoodchemical.com/ProductsList.aspx?{}'
enamine_url = 'https://www.enaminestore.com/api?code={}&currency=USD'

def get_sigma_aldrich_info(vid):
    res = url_request(
        sigma_aldrich_url.format(vid), 
        post = 'loadFor=PRD_RS',
        headers = {'Cookie': 'SialLocaleDef=WebLang~-1|CountryCode~CA|'},
    )

    # get prices
    cheapest_per_gram = 1e15
    prices = re.findall(r"<INPUT id='currencyAndPriceInPDP([^>]*)' name='currencyAndPriceInPDP([^>]*)' TYPE='hidden' Value=\"([^>]*)\">", res)
    for price_set in prices:
        _, qty_str = price_set[0].split('-')
        qty, unit = re.match(r'(\d+)([a-zA-Z]+)', qty_str).groups()
        qty = float(qty)
        if unit == 'MG':
            qty /= 1000
        currency, value = price_set[2].split('#@#')
        value = float(value.replace(',', ''))
        if value/qty < cheapest_per_gram:
            cheapest_per_gram = value/qty

    # get shipping
    date_fmt = '%m/%d/%y'
    earliest_date = datetime(2099, 12, 31, 0, 0)
    dates = re.findall(r"on (\d\d/\d\d/\d\d)", res)
    for d in dates:
        dtmp = datetime.strptime(d, date_fmt)
        if dtmp < earliest_date:
            earliest_date = dtmp

    return 0.76*cheapest_per_gram, earliest_date.strftime(date_fmt) 



def get_oakwood_info(vid):
    res = url_request( oakwood_url.format(vid), timeout=20)

    cheapest_per_gram = 1e15
    availability = 'unknown'

    # get prices
    prices = re.findall(r'<span id="_ctl0_ContentPlaceHolder1_MyGrid__ctl(.*)_([^<]*)">(.*)</span>', res)

    data = {}
    for price in prices:
        idx = price[0]
        if data.get(idx, None) is None:
            data[idx] = {'price': 0., 'qty': 0.}

        if price[1][:-1] == 'ProductID':
            qty, unit = re.match(r'(\d+)([a-zA-Z]+)', price[2].split('-')[1]).groups()
            qty = float(qty)
            if unit == 'mg':
                qty /= 1000
            data[idx]['qty'] = qty

        if price[1][:-1] == 'PriceBox':
            data[idx]['price'] = float(price[2][1:])

        if price[1] == 'InStock':
            availability = 'In Stock'

        if price[1] == 'OutOfStock' and availability == 'unknown':
            availability = 'Out of Stock'

    for d, v in data.items():
        if v['price']/v['qty'] < cheapest_per_gram:
            cheapest_per_gram = v['price']/v['qty']

    return cheapest_per_gram, availability

def get_enamine_info(vid):
    data = url_request( 
        enamine_url.format(vid), 
        headers={'Authorization': 'login=<your login info>'}, 
        fmt='JSON'
    )

    cheapest_per_gram = 1e15
    availability = 'Out of Stock'

    for entry in data['data']:
        if entry['packs'] is not None:
            for pack in entry['packs']:
                unit = pack['measure']
                qty = float(pack['amount'])
                price = float(pack['price'])
                if unit == 'mg':
                    qty /= 1e3
                if price/qty < cheapest_per_gram:
                    cheapest_per_gram = price/qty
    
            if float(entry['availability']) > 0.0:
                availability = 'In Stock'

    time.sleep(30)
    return cheapest_per_gram, availability


def get_sigma_aldrich_cheapest_old(vid):
    res = url_request(
        sigma_aldrich_url.format(vid), 
        post = 'loadFor=PRD_RS',
        headers = {'Cookie': 'SialLocaleDef=WebLang~-1|CountryCode~CA|'},
    )

    # get prices
    cheapest = 1e12
    prices = re.findall(r"<INPUT id='currencyAndPriceInPDP([^>]*)' name='currencyAndPriceInPDP([^>]*)' TYPE='hidden' Value=\"([^>]*)\">", res)
    for price_set in prices:
        qty_str = price_set[0].split('-')[1]
        qty, unit = re.match(r'(\d+)([a-zA-Z]+)', qty_str).groups()
        qty = float(qty)
        if unit == 'MG':
            qty /= 1000
        if unit == 'EA':
            return cheapest, 'CPR'

        currency, value = price_set[2].split('#@#')
        value = float(value.replace(',', ''))
        if value < cheapest:
            cheapest = value

    # get shipping
    date_fmt = '%m/%d/%y'
    earliest_date = datetime(2099, 12, 31, 0, 0)
    dates = re.findall(r"on (\d\d/\d\d/\d\d)", res)
    for d in dates:
        dtmp = datetime.strptime(d, date_fmt)
        if dtmp < earliest_date:
            earliest_date = dtmp

    return 0.76*cheapest, earliest_date.strftime(date_fmt), 


def get_sigma_aldrich_cheapest(vid):
    sa_url = "https://www.sigmaaldrich.com/CA/en/product/ALDRICH/{}"

    options = webdriver.firefox.options.Options()
    options.headless = True

    with webdriver.Firefox(options=options) as driver:
        driver.get(sa_url.format(vid))

        # find weight and id
        cheapest_per_gram = 1e15

        elements = driver.find_elements_by_css_selector("tr")
        for e in elements:
            pid = e.get_property('id').strip()
            if len(pid) <= 0:
                continue

            # get price
            pid_words = pid.split('-')
            price_id = '-'.join(pid_words[:2] + ['price'] + pid_words[2:])
            price_ele = driver.find_element_by_id(price_id)
            price = float(price_ele.text.split('$')[1]) * 0.79

            # parse weight
            qty_str = pid_words[-1]
            qty, unit = re.match(r'(\d+)([a-zA-Z]+)', qty_str).groups()

            qty = float(qty)
            if unit == 'MG':
                qty /= 1000

            # update cheapest price
            price_per_gram = price/qty
            if price_per_gram < cheapest_per_gram:
                cheapest_per_gram = price_per_gram

        date_fmt = '%m/%d/%y'
        availability = date.today().strftime(date_fmt)

    print(cheapest_per_gram, availability)
    return cheapest_per_gram, availability


def get_oakwood_cheapest_old(vid):
    res = url_request( oakwood_url.format(vid), timeout=20)

    cheapest = 1e8
    availability = 'unknown'

    # get prices
    prices = re.findall(r'<span id="_ctl0_ContentPlaceHolder1_MyGrid__ctl(.*)_([^<]*)">(.*)</span>', res.decode('utf-8'))

    data = {}
    for price in prices:
        idx = price[0]
        if data.get(idx, None) is None:
            data[idx] = {'price': 0., 'qty': 0.}

        if price[1][:-1] == 'ProductID':
            qty, unit = re.match(r'(\d+)([a-zA-Z]+)', price[2].split('-')[1]).groups()
            qty = float(qty)
            if unit == 'mg':
                qty /= 1000
            data[idx]['qty'] = qty

        if price[1][:-1] == 'PriceBox':
            data[idx]['price'] = float(price[2][1:])

        if price[1] == 'InStock':
            availability = 'In Stock'

        if price[1] == 'OutOfStock' and availability == 'unknown':
            availability = 'Out of Stock'

    for d, v in data.items():
        if v['price'] < cheapest:
            cheapest = v['price']

    return cheapest, availability

def get_oakwood_cheapest(vid):
    res = url_request( oakwood_url.format(vid), timeout=20)

    cheapest = 1e8
    availability = 'unknown'

    # get prices
    prices = re.findall(r'<span id="_ctl0_ContentPlaceHolder1_MyGrid__ctl(.*)_([^<]*)">(.*)</span>', res)


    data = {}
    for price in prices:
        idx = price[0]
        if data.get(idx, None) is None:
            data[idx] = {'price': 0., 'qty': 0.}

        if price[1][:-1] == 'ProductID':
            qty, unit = re.match(r'(\d+)([a-zA-Z]+)', price[2].split('-')[1]).groups()
            qty = float(qty)
            if unit == 'mg':
                qty /= 1000
            data[idx]['qty'] = qty

        if price[1][:-1] == 'PriceBox':
            data[idx]['price'] = float(price[2][1:])

        if price[1] == 'InStock':
            availability = 'In Stock'

        if price[1] == 'OutOfStock' and availability == 'unknown':
            availability = 'Out of Stock'

    for d, v in data.items():
        if v['price'] < cheapest:
            cheapest = v['price']

    print(cheapest, availability)
    return cheapest, availability




def get_vendor_and_vid(url):
    vendor, text = re.match(r'\w+\://[-a-zA-Z]+\.([-a-zA-Z]+)\.com(.+)', url).groups()
    vid = None
    if vendor == 'sigmaaldrich':
        match = re.match(r'/\w+/\w+/\w+/(.+)', text)
        vid = match.groups()[0].split('?')[0]
    elif vendor == 'oakwoodchemical':
        match = re.match(r'/[^?]+\?(.+)\&Ext', text)
        if match is None:
            vid = None
        else:
            vid = match.groups()[0]
    elif vendor == 'enaminestore':
        match = re.match(r'/catalog/(.+)', text)
        vid = match.groups()[0]

    return vendor, vid


vendor_parsefunc_dict = {
    'sigmaaldrich': get_sigma_aldrich_cheapest,
    'oakwoodchemical': get_oakwood_cheapest,
}

def get_price(long_url):
    urls = re.findall(r'http[^!$]*', long_url)

    data = {k: {'price': 1e15, 'availability': '', 'url': ''} for k in vendor_parsefunc_dict.keys()}
    for url in urls:
        vendor, vid = get_vendor_and_vid(url)
        print([vendor, vid])

        get_info = vendor_parsefunc_dict.get(vendor, None)

        if get_info is not None:
            try:
                price, availability = get_info(vid)
                if price < data[vendor]['price']:
                    data[vendor]['price'] = price
                    data[vendor]['availability'] = availability
                    data[vendor]['url'] = url
            except:
                print('TIMEOUT')
                pass

    return data

def get_vendor_prices(long_url):
    data = get_price(long_url)
    price_list = []
    for k, v in data.items():
        price_dict = {}
        if k == 'sigmaaldrich':
            price_dict['name'] = 'Sigma'
        elif k == 'oakwoodchemical':
            price_dict['name'] = 'Oakwood'
        elif k == 'enaminestore':
            price_dict['name'] = 'Enamine'
        else:
            continue
            
        price_dict['price'] = round(v['price'],2) if v['price'] < 1e7 else 'Nan'
        price_dict['availability'] = v['availability'] if len(v['availability']) > 0 else 'Nan'
        if price_dict['price'] == 'Nan':
            price_dict['url'] = ''
        else:
            price_dict['url'] = v['url']
        price_list.append(price_dict)
    return price_list

