#!/usr/bin/env python
# encoding: utf-8
"""
A settings file that is initialised by the user to configure the backend class.

Should be easily generalised to other settings that have similar requirements.
"""

def init(backend_obj):
    global be
    be = backend_obj