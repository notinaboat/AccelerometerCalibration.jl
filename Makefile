PACKAGE := $(shell basename $(PWD))
README.md: src/$(PACKAGE).jl
	julia --project -e "using $(PACKAGE); println(Docs.doc($(PACKAGE)))" > $@
