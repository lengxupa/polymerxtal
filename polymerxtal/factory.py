from .trajectoryadapter import TrajectoryAdapter
from .potential import LJ, Buckingham

_toolkits = {}

def register(toolkit_name, toolkit_class):
    if not issubclass(toolkit_class, TrajectoryAdapter):
        raise TypeError('{0} is not a TrajectoryAdapter'.format(toolkit_class))
    _toolkits[toolkit_name] = toolkit_class

def trajectory_factory(trajectory_toolkit, **kwargs):
    if trajectory_toolkit not in _toolkits.keys():
       raise TypeError('Toolkit not found')
       
    traj_analysis = _toolkits[trajectory_toolkit](kwargs['filename'])
    
    return traj_analysis

def potential_factory(potential_type, **kwargs):
   cls_dict = dict(LJ=LJ, Buckingham=Buckingham)

   if potential_type not in cls_dict.keys():
       raise Exception('Potential type not found')

   cls = cls_dict[potential_type]

   cls_instance = cls(**kwargs)
   return cls_instance
