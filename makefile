all: cython

cython:
	python2 setup.py build_ext --inplace
	$(MAKE) -C xlmhg
