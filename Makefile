all:
	python3 sgp4_op.py; \
	cd ADCS-Codes; \
	g++ -c adcs_engine_new.cpp; \
	g++ -c ephemerides_new.cpp; \
	g++ -c igrfnew.cpp; \
	g++ -c frame_conversion_new.cpp; \
	g++ -c OP_new.cpp; \
	g++ -c programs_new.cpp; \
	g++ -c sun_model_new.cpp; \
	g++ -c estimation.cpp; \
	g++ -c disturbance_torque.cpp; \
	g++ -c bdot.cpp; \
	g++ -c pid.cpp; \
	g++ -c integrator.cpp; \
	g++ -c reaction_wheel.cpp; \
	g++ -c power_gen.cpp; \
	g++ -c magnetorquer_torque.cpp; \
	g++ -o test adcs_engine_new.o ephemerides_new.o igrfnew.o frame_conversion_new.o OP_new.o programs_new.o sun_model_new.o estimation.o reaction_wheel.o magnetorquer_torque.o power_gen.o disturbance_torque.o bdot.o integrator.o pid.o; \
	./test; \
	cd ../py-igrf/release; python3 pyIGRF.py; \
    cd ../../plot_variables; python3 plot_variables.py; 

# -fsanitize=address -static-libasan
