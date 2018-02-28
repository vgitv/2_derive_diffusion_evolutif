# -----------------------------------------------------------------------------------------------------------
#  variables
# -----------------------------------------------------------------------------------------------------------
# nom projet
PROJECT = 2_derive_diffusion_evolutif

#
SOURCES = main.f95 math.f95
MODULES =          math.mod
OBJECTS = main.o   math.o
CC      = gfortran -fbounds-check
EXEC    = truc
OTHER   =
PDF     =



# -----------------------------------------------------------------------------------------------------------
#  compilation rules
# -----------------------------------------------------------------------------------------------------------
# exécutable
$(EXEC) : $(OBJECTS)
	$(CC) $(OBJECTS) -o $(EXEC)

# programme principal
main.o : main.f95 $(MODULES)
	$(CC) -c main.f95 -o main.o

# modules
math.o math.mod : math.f95
	$(CC) -c math.f95 -o math.o



# -----------------------------------------------------------------------------------------------------------
# phony targets
# -----------------------------------------------------------------------------------------------------------
# exécuter l'exécutable (utile pour utiliser F5 dans vim)
.PHONY : make
make :
	./$(EXEC)

# supprimer les fichiers objet et l'exécutable s'ils existent
.PHONY : clean
clean :
	rm -f $(EXEC) $(OBJECTS) $(MODULES)

# effacer le contenu des dossiers d'entrees et de sorties
.PHONY : del
del :
	rm -f sorties/* graphes/*

# ouvrir les fichiers du projet dans des onglets de vim
.PHONY : open
open :
	vim -p $(SOURCES) Makefile Plot.gnu $(OTHER)

# tout compiler et lancer gdb (segmentation fault)
.PHONY : gdb
gdb :
	$(CC) -g $(SOURCES) -o $(EXEC) -lm && gdb ./$(EXEC)

# clean et tarer le dossier
.PHONY : tar
tar :
	make clean
	rm -f .*.swp
	tar -zcvf ../$(PROJECT).tar.gz ../$(PROJECT)

# sauvegarder ancienne version
.PHONY : save
save :
	make clean
	rm -rf ../old_$(PROJECT)
	cp -r ../$(PROJECT) ../old_$(PROJECT)

.PHONY : pdf
pdf :
	xdg-open $(PDF)

#
.PHONY : clean
coffe :
	@echo "  (\n   )\n c[]"
