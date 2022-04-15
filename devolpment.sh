#!/bin/bash
mkdir -p static/svg

if [[ -z "$(~/a2g2-singularity/a2g2 pg list | grep a2g2_db)" ]]; then
    ~/a2g2-singularity/a2g2 pg start
fi
~/a2g2-singularity/a2g2 exec "unset DJANGO_SETTINGS_MODULE; ./manage.py runserver 7000"
