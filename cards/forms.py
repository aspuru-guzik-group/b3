from django import forms

class MolForm(forms.Form):
    smiles = forms.CharField(
        max_length=200, 
        widget=forms.TextInput(attrs={'class': 'form-control'})
    )
    info = forms.CharField(
        max_length=500, 
        required = False,
        widget=forms.Textarea(attrs={'class': 'form-control', 'rows': 3})
    )
    prop = forms.CharField(
        max_length=500, 
        required = False,
        widget=forms.Textarea(attrs={'class': 'form-control', 'rows': 3})
    )

    def disable_smiles(self):
        self.fields['smiles'].widget.attrs['readonly'] = True

class BlockForm(forms.Form):
    smiles = forms.CharField(
        max_length=1500, 
        required = False,
        widget=forms.TextInput(attrs={'class': 'form-control'})
    )

    smiles_list = forms.CharField(
        required = False,
        widget=forms.Textarea(attrs={'class': 'form-control', 'rows': 3})
    )

    cas_list = forms.CharField(
        required = False,
        widget=forms.Textarea(attrs={'class': 'form-control', 'rows': 3})
    )

    button_type = forms.CharField(
        required=True,
        initial='add',
        widget=forms.HiddenInput(),
    )
    
