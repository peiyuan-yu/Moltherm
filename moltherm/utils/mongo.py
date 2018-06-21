from pymongo import MongoClient


def mongo_collection(host, port, db_name, collection_name,
                     username=None, password=None):
    """
    Get the collection of a dataset

    """
    if username and password:
        mongo_uri = 'mongodb://{}:{}@{}:{}/{}'.format(username, password,
                                                      host, port, db_name)
        client = MongoClient(mongo_uri)
    else:
        client = MongoClient(host, port)
    db = client[db_name]
    collection = db[collection_name]
    return collection
