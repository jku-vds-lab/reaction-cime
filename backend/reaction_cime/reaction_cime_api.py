import pickle
from flask import Blueprint, request
import logging

_log = logging.getLogger(__name__)

reaction_cime_api = Blueprint('reaction_cime', __name__)

@reaction_cime_api.route('/hello', methods=['GET'])
def get_uploaded_files_list():
    return "Hello World"
