package org.sbpo2025.challenge;

import ilog.concert.IloException;
import ilog.concert.IloIntVar;
import ilog.concert.IloLinearIntExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;
import org.apache.commons.lang3.time.StopWatch;

import java.util.*;
import java.util.concurrent.TimeUnit;

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
    protected IloCplex prob;

    public ChallengeSolver(
            List<Map<Integer, Integer>> orders, List<Map<Integer, Integer>> aisles, int nItems, int waveSizeLB, int waveSizeUB) {
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
        double valorObjetivoActual;
        double rango_k = Math.log(this.waveSizeUB - ((double) this.waveSizeLB / cantPasillos)) + Math.log(1 / EPSILON);
        Map<Integer, Boolean> dictWSol;
        Map<Integer, Boolean> dictASol;

        if (cantPasillos <= Math.ceil(rango_k)) {
            double maximo = INF;
            double valorObjetivo;
            Map<Integer, Boolean>[] solucionActual;
            for (int aPrima = 1;  aPrima < cantPasillos + 1 ; aPrima++){ //incluyo ese máximo
                solucionActual = planteoPasillosFijos(aPrima);
                resolver_lp(prob);
                valorObjetivoActual = prob.getObjValue() / aPrima;
                if (maximo < valorObjetivoActual){
                    maximo = valorObjetivoActual;
                    dictWSol = solucionActual[0];
                    dictASol = solucionActual[1];
                    valorObjetivo = valorObjetivoActual;
                    System.out.println(valorObjetivo);
                }
            }
        } else {
            int l = this.waveSizeLB;
            int u = this.waveSizeUB;
            int k;
            Map<Integer, Boolean>[] dictSol = new Map[2];
            while (u - l > EPSILON) {
                k = Math.floorDiv(u - l, 2);
                dictSol = planteo_busqueda_binaria(k);
                resolver_lp(prob);
                valorObjetivoActual = prob.getObjValue();
                if (valorObjetivoActual > 0) { //Era infactible con ese valor de k alto
                    u = k;
                }else{
                    l = k;
                }
            }
            dictWSol = dictSol[0];
            dictASol = dictSol[1];
            System.out.println(prob.getObjValue());
        }

        ChallengeSolution solution = new ChallengeSolution(dictWSol.keySet().stream().filter(k -> dictWSol.get(k)), dictASol.keySet().stream().filter(k -> dictASol.get(k));
        return null;
    }

    private List<IloNumVar> planteoPasillosFijos(int aPrima) throws IloException {
        List<IloNumVar> resPasillos = new ArrayList<>();
        IloNumVar[] listaW1 = new IloNumVar[this.orders.size()];
        IloNumVar[] listaA1 = new IloNumVar[this.aisles.size()];

        //Si la orden pertenece a la wave
        for(int o = 0; o < this.orders.size(); o++){
            listaW1[o] = prob.boolVar(String.format("W_%d",o));
        }

        //Si se usa el pasillo a
        for(int a = 0; a < aPrima; a++){
            listaA1[a] = prob.boolVar(String.format("A_%d",a));
        }

        //La cantidad de elementos que se toman en la Wave está en rango
        IloLinearIntExpr suma = prob.linearIntExpr();
        for (int o = 0; o < this.orders.size(); o++) {
            for (Map.Entry<Integer, Integer> i : this.orders.get(o).entrySet()) {
                suma.addTerm(i.getValue().intValue(), (IloIntVar) listaW1[o]);
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
                izq.addTerm(u_oi, (IloIntVar) listaW1[o]);
            }
            for (int a = 0; a < this.aisles.size(); a++) {
                int u_ai = this.aisles.get(a).getOrDefault(i, 0);
                der.addTerm(u_ai, (IloIntVar) listaA1[a]);
            }
            prob.addLe(izq, der); // Agrega la restricción
        }

        //La cantidad de pasillos usados es A* (pasado por parámetro)
        IloLinearIntExpr sumaDeA = prob.linearIntExpr();
        for (int a = 0; a < this.aisles.size(); a++) {
            sumaDeA.addTerm(1, (IloIntVar) listaA1[a]);
        }
        prob.addEq(sumaDeA, aPrima);
        prob.addMaximize(suma);

        Collections.addAll(resPasillos, listaW1);
        Collections.addAll(resPasillos, listaA1);
        return resPasillos;
    }

    private List<IloNumVar> planteo_busqueda_binaria(int k) throws IloException {
        List<IloNumVar> resBinaria = new ArrayList<>();
        IloNumVar[] listaW2 = new IloNumVar[this.orders.size()];
        IloNumVar[] listaA2 = new IloNumVar[this.aisles.size()];

        //Si la orden pertenece a la wave
        for(int o = 0; o < this.orders.size(); o++){
            listaW2[o] = prob.boolVar(String.format("W_%d",o));
        }

        //Si se usa el pasillo a
        for(int a = 0; a < this.aisles.size(); a++){
            listaA2[a] = prob.boolVar(String.format("A_%d",a));
        }

        //La cantidad de elementos que se toman en la Wave está en rango
        IloLinearIntExpr suma = prob.linearIntExpr();
        for (int o = 0; o < this.orders.size(); o++) {
            for (Map.Entry<Integer, Integer> i : this.orders.get(o).entrySet()) {
                suma.addTerm(i.getValue().intValue(), (IloIntVar) listaW2[o]);

        //Restricción de la forma f = k * g
        IloLinearIntExpr pasillosUtilizados = prob.linearIntExpr();
            for (int a = 0; a < this.aisles.size(); a++) {
                pasillosUtilizados.addTerm(k, (IloIntVar) listaA2[a]);
            }

        prob.addLe(pasillosUtilizados.addTerm(-1, EPSILON), suma);
        prob.addGe(pasillosUtilizados.addTerm(1, EPSILON), suma);

        //Función objetivo: REVISAR
        for (int i = 0; i < this.nItems; i++) {
            IloLinearIntExpr izq = prob.linearIntExpr();
            IloLinearIntExpr der = prob.linearIntExpr();

            for (int o = 0; o < this.orders.size(); o++) {
                int u_oi = this.orders.get(o).getOrDefault(i, 0);
                izq.addTerm(u_oi, (IloIntVar) listaW2[o]);
            }
            for (int a = 0; a < this.aisles.size(); a++) {
                int u_ai = this.aisles.get(a).getOrDefault(i, 0);
                der.addTerm(u_ai, (IloIntVar) listaA2[a]);
            }
            prob.addMinimize(izq.addTerm(-1, der), "minimize");
        }

        Collections.addAll(resBinaria, listaW2);
        Collections.addAll(resBinaria, listaA2);
        return resBinaria;
    }

    private void resolver_lp(){
        //Definir los parametros del solver
        prob.setParam(IloCplex.Param.MIP.Tolerances.AbsMIPGap, TOLERANCE);
        //Resolver el lp
        prob.solve();
    }
    /*
     * Get the remaining time in seconds
     */
    protected long getRemainingTime(StopWatch stopWatch) {
        return Math.max(
                TimeUnit.SECONDS.convert(MAX_RUNTIME - stopWatch.getTime(TimeUnit.MILLISECONDS), TimeUnit.MILLISECONDS),
                0);
    }

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
}
