# -*- Makefile -*-
# ----------------------------------------------------------------------

TARGETS = \
  $(SEQDB_LIB) \
  $(SEQDB_PY_LIB) \
  $(DIST)/seqdb3-scan \
  $(DIST)/seqdb3 \
  $(DIST)/seqdb3-chart-sequenced \
  $(DIST)/seqdb3-compare-sequences \
  $(DIST)/seqdb3-chart-clades \
  $(DIST)/seqdb3-seqid-by-name \
  $(DIST)/seqdb3-chart-dates \
  $(DIST)/test-insertions-deletions

SEQDB_SOURCES = \
  seqdb.cc \
  scan-fasta.cc ncbi.cc scan-sequence.cc scan-align.cc eliminate-identical.cc scan-deletions.cc scan-lineages.cc scan-match-hidb.cc create.cc \
  seq-id.cc seqdb-parse.cc compare.cc aa-at-pos.cc

SEQDB_LIB_MAJOR = 3
SEQDB_LIB_MINOR = 0
SEQDB_LIB_NAME = libseqdb
SEQDB_LIB = $(DIST)/$(call shared_lib_name,$(SEQDB_LIB_NAME),$(SEQDB_LIB_MAJOR),$(SEQDB_LIB_MINOR))

# ----------------------------------------------------------------------

all: install

include $(ACMACSD_ROOT)/share/Makefile.config

LDLIBS = \
  $(AD_LIB)/$(call shared_lib_name,libacmacsbase,1,0) \
  $(AD_LIB)/$(call shared_lib_name,liblocationdb,1,0) \
  $(AD_LIB)/$(call shared_lib_name,libacmacsvirus,1,0) \
  $(AD_LIB)/$(call shared_lib_name,libacmacschart,2,0) \
  $(AD_LIB)/$(call shared_lib_name,libhidb,5,0) \
  $(XZ_LIBS) $(PYTHON_LIBS) $(CXX_LIBS)

# $(AD_LIB)/$(call shared_lib_name,libacmacschart,2,0) \

# ----------------------------------------------------------------------

install: install-headers $(TARGETS)
	$(call install_lib,$(SEQDB_LIB))
	$(call symbolic_link_wildcard,$(DIST)/seqdb*,$(AD_BIN))

test: install
	test/test
.PHONY: test

# ----------------------------------------------------------------------

$(SEQDB_LIB): $(patsubst %.cc,$(BUILD)/%.o,$(SEQDB_SOURCES)) | $(DIST) install-headers
	$(call echo_shared_lib,$@)
	$(call make_shared_lib,$(SEQDB_LIB_NAME),$(SEQDB_LIB_MAJOR),$(SEQDB_LIB_MINOR)) $(LDFLAGS) -o $@ $^ $(LDLIBS)

$(DIST)/%: $(BUILD)/%.o | $(SEQDB_LIB)
	$(call echo_link_exe,$@)
	$(CXX) $(LDFLAGS) -o $@ $^ $(SEQDB_LIB) $(LDLIBS) $(AD_RPATH)

# ======================================================================
### Local Variables:
### eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
### End:
