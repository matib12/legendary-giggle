FROM gcc:12-bullseye

WORKDIR /app

COPY functions.F90 nbsimulation.F90 ./

CMD sh -c "gfortran -c functions.F90 && \
           gfortran functions.o nbsimulation.F90 -o simulation && \
           rm *.o *.mod && \
           mkdir -p results && \
           mv input.TXT results/ 2>/dev/null || cp input.TXT results/ && \
           cd results && \
           ../simulation"
