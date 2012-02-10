## Config
CXX = g++
CFLAGS = -g #-O2
INCLUDE = /usr/include # $(CUBA)
CPPFLAGS = -I$(HOME)/local/include -I$(INCLUDE) -I. 
LDFLAGS = -g -L$(HOME)/local/lib
SHELL = /bin/sh


## Files


APP= main

SRCS = bg gating spannung solve my_commands \
	$(addsuffix .cc, $(basename $(APP) ) )

## Rules
all: $(APP)

OBJS = $(addsuffix .o, $(basename $(SRCS)))
DEPS = $(addprefix .deps/, $(addsuffix .P, $(basename $(SRCS))))

#$(APP): .deps $(OBJS) $(GINAC_BUILD)/ginac/libginac.la
#	$(SHELL) libtool --mode=link $(CXX) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

$(APP): .deps $(OBJS)
	$(CXX) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)
#	$(SHELL) /usr/bin/libtool --mode=link $(CXX) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

%.o: %.cc
	$(CXX)  -I$(INCLUDE) $(CFLAGS) $(CPPFLAGS) -Wp,-MD,.deps/$(*F).P -o $@ -c $<

clean:
	rm -f $(OBJS) $(APP)
	rm -rf .libs *.o

distclean: clean
	rm -rf .deps *.o *~
.deps:
	@[ -d .deps ] || mkdir .deps >/dev/null 2>&1

-include $(DEPS)


# Makefile adapted by Markus Knodel 2009 2010 2011


