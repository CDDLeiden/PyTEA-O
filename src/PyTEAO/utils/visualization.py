import importlib
import pkgutil
import PyTEAO.visualization.subplots

SUBPLOT_REGISTRY = {}

def register(name):
	"""Decorator to register a subplot component class."""
	def decorator(cls):
		if name in SUBPLOT_REGISTRY:
			raise ValueError(f"Duplicate subplot key: '{name}'")
		SUBPLOT_REGISTRY[name] = cls
		return cls
	return decorator

def import_all_subplots():
	"""Dynamically import all subplot modules to populate the registry."""
	package = PyTEAO.visualization.subplots
	for _, module_name, _ in pkgutil.iter_modules(package.__path__):
		importlib.import_module(f"{package.__name__}.{module_name}")
