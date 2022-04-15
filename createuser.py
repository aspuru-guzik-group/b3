import os, django
os.environ['DJANGO_SETTINGS_MODULE'] = 'bbb.settings'
django.setup()

from django.contrib.auth.models import User
user = User.objects.create_user('name', password='password')
user.save()
