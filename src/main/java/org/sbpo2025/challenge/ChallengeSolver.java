package org.sbpo2025.challenge;

import java.util.*;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import ilog.concert.*;
import org.apache.commons.lang3.time.StopWatch;

import ilog.cplex.IloCplex;

public class ChallengeSolver {
    private final long MAX_RUNTIME = 600000; // milliseconds; 10 minutes
    private final double EPSILON = 0.001;
    private final double INF = Double.MIN_VALUE;
    private final double TOLERANCE = Math.exp(-6);

    protected List<Map<Integer, Integer>> orders;
    protected List<Map<Integer, Integer>> aisles;
    protected int nItems;
    protected int waveSizeLB;
    protected int waveSizeUB;

    public ChallengeSolver(
            List<Map<Integer, Integer>> orders, List<Map<Integer, Integer>> aisles, int nItems, int waveSizeLB, int waveSizeUB) throws IloException {
        this.orders = orders;
        this.aisles = aisles;
        this.nItems = nItems;
        this.waveSizeLB = waveSizeLB;
        this.waveSizeUB = waveSizeUB;
    }

    public ChallengeSolution solve(StopWatch stopWatch) throws IloException {
        // DEFINICIÓN DEL MODELO Y SOLVER
        //Decidimos qué modelo vamos a usar (el que tenga que resolver menos PL
        int cantPasillos = this.aisles.size();
        double epsilon = 1/ (double) cantPasillos;
        double rango_k = (Math.log(this.waveSizeUB - ((double) this.waveSizeLB / cantPasillos)) - Math.log(epsilon))/(Math.log(2.0));
        List<Boolean> dictWSol = List.of();
        List<Boolean> dictASol = List.of();
        IloCplex prob = new IloCplex();
        prob.setOut(null);

        System.out.println(String.format("Cantidad de pasillos: %d", cantPasillos));
        System.out.println(String.format("LB: %d", this.waveSizeLB));
        System.out.println(String.format("UB: %d", this.waveSizeUB));

        System.out.println(this.waveSizeUB - ((double) this.waveSizeLB / cantPasillos));
        System.out.println(Math.log(this.waveSizeUB - ((double) this.waveSizeLB / cantPasillos)));
        System.out.println(rango_k);
        // TODO: Debuggear binaria
//        List<List<Boolean>> solucionActual = planteo_busqueda_binaria(prob, epsilon);
//            dictWSol = solucionActual.get(0);
//            dictASol = solucionActual.get(1);
//            System.out.println(prob.getObjValue());

        if (true || cantPasillos <= rango_k) {
            List<List<Boolean>> solucionActual = planteoPasillosFijos(prob);
            dictWSol = solucionActual.get(0);
            dictASol = solucionActual.get(1);
        }/*
        else {
            System.out.println("Eligio binaria");
            List<List<Boolean>> solucionActual = planteo_busqueda_binaria(prob, epsilon);
            dictWSol = solucionActual.get(0);
            dictASol = solucionActual.get(1);
            System.out.println(prob.getObjValue());
        }*/

        List<Boolean> finalDictWSol = dictWSol;
        List<Boolean> finalDictASol = dictASol;
        Set<Integer> finalOrder = IntStream.range(0, orders.size()).filter(finalDictWSol::get).boxed().collect(Collectors.toSet());
        Set<Integer> finalAisle = IntStream.range(0, aisles.size()).filter(finalDictASol::get).boxed().collect(Collectors.toSet());

        return new ChallengeSolution(finalOrder, finalAisle);
    }

    private List<List<Boolean>> planteoPasillosFijos(IloCplex prob) throws IloException {
        List<List<Boolean>> resPasillos = new ArrayList<>();

        if (this.orders.isEmpty() || this.aisles.isEmpty()) {
            throw new IloException("Error: No hay órdenes o pasillos disponibles.");
        }

        IloIntVar[] listaW1 = new IloIntVar[this.orders.size()];
        IloIntVar[] listaA1 = new IloIntVar[this.aisles.size()];

        //Si la orden pertenece a la wave
        for (int o = 0; o < this.orders.size(); o++) {
            listaW1[o] = prob.boolVar(String.format("W_%d", o));
        }

        //Si se usa el pasillo a
        for (int a = 0; a < this.aisles.size(); a++) {
            listaA1[a] = prob.boolVar(String.format("A_%d", a));
        }

        //La cantidad de elementos que se toman en la Wave está en rango
        IloLinearIntExpr suma = prob.linearIntExpr();
        for (int o = 0; o < this.orders.size(); o++) {
            for (Map.Entry<Integer, Integer> i : this.orders.get(o).entrySet()) {
                suma.addTerm(i.getValue(), listaW1[o]);
            }
        }
        prob.addLe(suma, this.waveSizeUB);
        prob.addGe(suma, this.waveSizeLB);

        // Todos los elementos que tomo de un pasillo tienen stock
        for (int i = 0; i < this.nItems; i++) {
            IloLinearIntExpr izq = prob.linearIntExpr();
            IloLinearIntExpr der = prob.linearIntExpr();

            for (int o = 0; o < this.orders.size(); o++) {
                int u_oi = this.orders.get(o).getOrDefault(i, 0);
                izq.addTerm(u_oi, listaW1[o]);
            }
            for (int a = 0; a < this.aisles.size(); a++) {
                int u_ai = this.aisles.get(a).getOrDefault(i, 0);
                der.addTerm(u_ai, listaA1[a]);
            }
            prob.addLe(izq, der); // Agrega la restricción
        }

        double maximo = INF;
        double valorObjetivoActual;
        List<Boolean> valoresW = List.of();
        List<Boolean> valoresA = List.of();
        prob.addMaximize(suma);
        int aPrima = 1;
        prob.setParam(IloCplex.Param.MIP.Tolerances.AbsMIPGap, TOLERANCE);
        IloLinearIntExpr sumaDeA = prob.linearIntExpr();
        for (int a = 0; a < this.aisles.size(); a++) {
            sumaDeA.addTerm(1, listaA1[a]);
        }
        while (aPrima < this.aisles.size() + 1 && maximo*aPrima <=  this.waveSizeUB) {
            //int aPrima = 2;
            //La cantidad de pasillos usados es A* (pasado por parámetro)

            IloConstraint restriccionA = prob.addEq(sumaDeA, aPrima);

            if (prob.solve()){
                //Resolver el lp
                valorObjetivoActual = prob.getObjValue() / aPrima;

                if (maximo <= valorObjetivoActual) {
                    maximo = valorObjetivoActual;

                    valoresW = new ArrayList<>();
                    valoresA = new ArrayList<>();

                    for (IloIntVar w : listaW1) {
                        valoresW.add(prob.getValue(w) >= 0.5);
                    }

                    for (IloIntVar a : listaA1) {
                        valoresA.add(prob.getValue(a) >= 0.5);
                    }
                }
                System.out.println(valorObjetivoActual);
                //prob.exportModel(String.format("model_%d.lp", aPrima));
            } else {
                System.out.println(String.format("Infactible para a'=%d", aPrima));
            }
            prob.remove(restriccionA);
            aPrima++;
        }

        resPasillos.add(valoresW);
        resPasillos.add(valoresA);

        return resPasillos;
    }


    private List<List<Boolean>> planteo_busqueda_binaria(IloCplex prob, double epsilon) throws IloException {
        List<List<Boolean>> resPasillos = new ArrayList<>();

        if (this.orders.isEmpty() || this.aisles.isEmpty()) {
            throw new IloException("Error: No hay órdenes o pasillos disponibles.");
        }

        IloIntVar[] listaW1 = new IloIntVar[this.orders.size()];
        IloIntVar[] listaA1 = new IloIntVar[this.aisles.size()];

        //Si la orden pertenece a la wave
        for (int o = 0; o < this.orders.size(); o++) {
            listaW1[o] = prob.boolVar(String.format("W_%d", o));
        }

        //Si se usa el pasillo a
        for (int a = 0; a < this.aisles.size(); a++) {
            listaA1[a] = prob.boolVar(String.format("A_%d", a));
        }

        IloNumVar z = prob.numVar(0, waveSizeUB, "z");

        //La cantidad de elementos que se toman en la Wave está en rango
        IloLinearIntExpr suma = prob.linearIntExpr();
        for (int o = 0; o < this.orders.size(); o++) {
            for (Map.Entry<Integer, Integer> i : this.orders.get(o).entrySet()) {
                suma.addTerm(i.getValue(), listaW1[o]);
            }
        }
        prob.addLe(suma, this.waveSizeUB);
        prob.addGe(suma, this.waveSizeLB);

        // Todos los elementos que tomo de un pasillo tienen stock
        for (int i = 0; i < this.nItems; i++) {
            IloLinearNumExpr izq = prob.linearNumExpr();
            IloLinearNumExpr der = prob.linearNumExpr();

            for (int o = 0; o < this.orders.size(); o++) {
                int u_oi = this.orders.get(o).getOrDefault(i, 0);
                izq.addTerm(u_oi, listaW1[o]);
            }
            for (int a = 0; a < this.aisles.size(); a++) {
                int u_ai = this.aisles.get(a).getOrDefault(i, 0);
                der.addTerm(u_ai, listaA1[a]);
            }
            der.addTerm(1, z);
            prob.addLe(izq, der); // Agrega la restricción
        }

        double maximo = INF;
        double valorObjetivoActual;
        List<Boolean> valoresW = List.of();
        List<Boolean> valoresA = List.of();
        prob.addMinimize(z);
        prob.setParam(IloCplex.Param.MIP.Tolerances.AbsMIPGap, TOLERANCE);
        IloLinearNumExpr sumaDeA = prob.linearNumExpr();
        for (int a = 0; a < this.aisles.size(); a++) {
            sumaDeA.addTerm(1, listaA1[a]);
        }

        double searchMin = (double) waveSizeLB / (double) this.aisles.size();
        double searchMax = waveSizeUB;
        double last_feasible_k = searchMax;

        while (searchMin + epsilon < searchMax){
            double k = (searchMax + searchMin) / 2;
            System.out.println(k);

            IloConstraint restriccionA = prob.addEq(prob.prod(k, sumaDeA), suma); // Convertirlo a restricciones ensanguchadas con epsilon

            boolean isSolved = prob.solve();

            if (isSolved){
                //Resolver el lp
                double z_obj = prob.getObjValue();
                last_feasible_k = k;
                System.out.println(String.format("Objetivo: %f", z_obj));

                if (z_obj > 0){
                    searchMax = k;
                }
                else{
                    searchMin = k;
                }
            } else {
                System.out.println(String.format("Infactible para k=%f", k));
                searchMin = k;
            }

            prob.remove(restriccionA);
        }
        
        prob.addEq(prob.prod(last_feasible_k, sumaDeA), suma);
        prob.solve();
        prob.exportModel("model_binaria.lp");

        valoresW = new ArrayList<>();
        valoresA = new ArrayList<>();

        for (IloIntVar w : listaW1) {
            valoresW.add(prob.getValue(w) >= 0.5);
        }

        for (IloIntVar a : listaA1) {
            valoresA.add(prob.getValue(a) >= 0.5);
        }

        resPasillos.add(valoresW);
        resPasillos.add(valoresA);

        return resPasillos;
    }
    /*
     * Get the remaining time in seconds
     */
    protected long getRemainingTime(StopWatch stopWatch) {
        return Math.max(
                TimeUnit.SECONDS.convert(MAX_RUNTIME - stopWatch.getTime(TimeUnit.MILLISECONDS), TimeUnit.MILLISECONDS),
                0);
    }
}
/*
    protected boolean isSolutionFeasible(ChallengeSolution challengeSolution) {
        Set<Integer> selectedOrders = challengeSolution.orders();
        Set<Integer> visitedAisles = challengeSolution.aisles();
        if (selectedOrders == null || visitedAisles == null || selectedOrders.isEmpty() || visitedAisles.isEmpty()) {
            return false;
        }

        int[] totalUnitsPicked = new int[nItems];
        int[] totalUnitsAvailable = new int[nItems];

        // Calculate total units picked
        for (int order : selectedOrders) {
            for (Map.Entry<Integer, Integer> entry : orders.get(order).entrySet()) {
                totalUnitsPicked[entry.getKey()] += entry.getValue();
            }
        }

        // Calculate total units available
        for (int aisle : visitedAisles) {
            for (Map.Entry<Integer, Integer> entry : aisles.get(aisle).entrySet()) {
                totalUnitsAvailable[entry.getKey()] += entry.getValue();
            }
        }

        // Check if the total units picked are within bounds
        int totalUnits = Arrays.stream(totalUnitsPicked).sum();
        if (totalUnits < waveSizeLB || totalUnits > waveSizeUB) {
            return false;
        }

        // Check if the units picked do not exceed the units available
        for (int i = 0; i < nItems; i++) {
            if (totalUnitsPicked[i] > totalUnitsAvailable[i]) {
                return false;
            }
        }

        return true;
    }

    protected double computeObjectiveFunction(ChallengeSolution challengeSolution) {
        Set<Integer> selectedOrders = challengeSolution.orders();
        Set<Integer> visitedAisles = challengeSolution.aisles();
        if (selectedOrders == null || visitedAisles == null || selectedOrders.isEmpty() || visitedAisles.isEmpty()) {
            return 0.0;
        }
        int totalUnitsPicked = 0;

        // Calculate total units picked
        for (int order : selectedOrders) {
            totalUnitsPicked += orders.get(order).values().stream()
                    .mapToInt(Integer::intValue)
                    .sum();
        }

        // Calculate the number of visited aisles
        int numVisitedAisles = visitedAisles.size();

        // Objective function: total units picked / number of visited aisles
        return (double) totalUnitsPicked / numVisitedAisles;
    }
}*/
//java -Djava.library.path=""HOME/Users/DAFNE/OneDrive/Documentos/GitHub/challenge-sbpo-2025\cplex\lib\cplex.jar"" -jar target/ChallengeSBPO2025-1.0.jar datasets/a/instance_0001.txt primerRes.txt
// java -Djava.library.path="C:/Users/DAFNE/challenge-sbpo-2025/cplex/bin/x64_win64" -jar target/ChallengeSBPO2025-1.0.jar datasets\a\instance_0001.txt salida.txt
// java -jar target/ChallengeSBPO2025-1.0.jar input_prueba.txt ./outputs.txt
