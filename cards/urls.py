from django.urls import path
from . import views

urlpatterns = [
    path('project', views.projects, name='projects'),
    path('project/<project>', views.projectpage, name='projects'),
    path('project/<project>/data.csv', views.projectpage_download, name='projects'),
    path('project/<project>/new', views.projectpage_addpage, name='projects'),
]
