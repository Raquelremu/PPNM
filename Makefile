MCS = mcs
MONO = mono

LIBRARY = sfuns.dll
MAIN_EXE = main.exe
OUTPUT_FILE = Out.txt

LIB_SRC = sfuns.cs
MAIN_SRC = main.cs

all: $(MAIN_EXE)

$(LIBRARY): $(LIB_SRC)
	$(MCS) -target:library -out:$(LIBRARY) $(LIB_SRC)

$(MAIN_EXE): $(MAIN_SRC) $(LIBRARY)
	$(MCS) -target:exe -out:$(MAIN_EXE) -reference:$(LIBRARY) $(MAIN_SRC)

run: $(MAIN_EXE)
	$(MONO) $(MAIN_EXE) > $(OUTPUT_FILE)
	cat $(OUTPUT_FILE)

clean:
	rm -f $(LIBRARY) $(MAIN_EXE) $(OUTPUT_FILE)

