# Generate dual pairs of algebras corresponding to Galois
# representations with image GL_2(F_3) or GL_2(F_5).

# DATAFILE must be set on the command line.  Use with results
# downloaded in text format from
# * l = 3:
# https://www.lmfdb.org/NumberField/?galois_group=8T23&hidecol=disc.galois.class_group_desc
# * l = 5:
# https://www.lmfdb.org/NumberField/?galois_group=24T1353&hidecol=disc.galois.class_group_desc
DATAFILE =

GP = gp -q -D parisizemax=4G

# Evaluate $(1) in GP with arguments $(2), redirecting standard output
# to $(3).
gp_eval_raw = echo $(1) | $(GP) $(2) > $(3)
gp_catch = "iferr("$(1)",e,write(\"/dev/stderr\",e);quit(1))"
# Evaluate $(1) in GP with arguments $(2), catching errors and writing
# the result to $(3).
gp_eval = $(call gp_eval_raw,"print("$(call gp_catch,$(1))");",$(2),$(3))

pairs:
	$(MAKE) $(shell sed -ne 's/^"\(.*\)".*/\1.gp/p' $(DATAFILE))

%.gp:
	$(call gp_eval,"to_lmfdb_format(dual_pair_from_GL2_field(Polrev("$(shell sed -ne "s/\"$*\".\(.*\)/\1/p" $(DATAFILE))")))",GL2.gp,$@)
