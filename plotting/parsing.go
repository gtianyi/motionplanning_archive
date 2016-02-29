package main

import (
	"os"
	"strings"
	"bufio"
	"fmt"
	"strconv"
)

type dataPoint struct {
	Solved         bool
	Length         float64
	Time           float64
	Precomputation float64
}

func parseParams(data []string, keys []string) map[string]string {
	params := map[string]string{}
	for _, line := range data {
		for _, key := range keys {
			vals := strings.Split(line, " = ")
			if vals[0] == key {
				params[key] = strings.TrimSpace(vals[1])
			}
		}
	}
	return params
}

func parseFBiasedRRTParams(data []string) map[string]string {
	keys := []string{"num_prm_edges", "omega", "prm_size", "state_radius", "sampler_initialization_time"}
	return parseParams(data, keys)
}

func parseFBiasedRRTShellParams(data []string) map[string]string {
	keys := []string{"num_prm_edges", "omega", "prm_size", "state_radius", "shell_preference", "shell_depth", "sampler_initialization_time"}
	return parseParams(data, keys)
}

func parsePlakuRRTParams(data []string) map[string]string {
	keys := []string{"num_prm_edges", "omega", "prm_size", "state_radius", "alpha", "b", "sampler_initialization_time"}
	return parseParams(data, keys)
}

func parseAstarParams(data []string) map[string]string {
	keys := []string{"num_prm_edges", "omega", "prm_size", "state_radius", "weight", "sampler_initialization_time"}
	return parseParams(data, keys)
}

func parseDijkstraParams(data []string) map[string]string {
	keys := []string{"num_prm_edges", "omega", "prm_size", "state_radius", "sampler_initialization_time"}
	return parseParams(data, keys)
}

func parseGreedyParams(data []string) map[string]string {
	keys := []string{"num_prm_edges", "omega", "prm_size", "state_radius", "sampler_initialization_time"}
	return parseParams(data, keys)
}

func parseSpeedyParams(data []string) map[string]string {
	keys := []string{"num_prm_edges", "omega", "prm_size", "state_radius", "use_d_or_e", "sampler_initialization_time"}
	return parseParams(data, keys)
}

func parseNewPlannerParams(data []string) map[string]string {
	keys := []string{"num_prm_edges", "prm_size", "state_radius", "sampler_initialization_time"}
	return parseParams(data, keys)
}

func nonRDBReader(filename string) (map[string]string, map[string][][]string, bool) {
	openFile, err := os.Open(filename)
	if err != nil {
		panic(err)
	}

	defer openFile.Close()

	scanner := bufio.NewScanner(openFile)
	data := []string{}
	capture := false
	for scanner.Scan() {
		str := scanner.Text()

		if capture {
			data = append(data, str)
		}

		capture = capture || (str == "1 planners")
		if str == "." {
			break
		}
	}

	params := map[string]string{}

	plannerName := strings.Replace(data[0], "control_", "", -1)

	if plannerName == "FBiased RRT" {
		params = parseFBiasedRRTParams(data)
	} else if plannerName == "FBiased RRT Shell" {
		params = parseFBiasedRRTShellParams(data)
	} else if plannerName == "Plaku RRT" {
		params = parsePlakuRRTParams(data)
	} else if plannerName == "Dijkstra" {
		params = parseDijkstraParams(data)
	} else if plannerName == "A*" {
		params = parseAstarParams(data)
		if params["weight"] != "1" {
			weight, _ := strconv.ParseFloat(params["weight"], 64)
			plannerName = fmt.Sprintf("WA* (w=%.2f)", weight)
		}
	} else if plannerName == "Speedy" {
		params = parseSpeedyParams(data)
		plannerName = fmt.Sprintf("Speedy (%s)", params["use_d_or_e"])
	} else if plannerName == "Greedy" {
		params = parseGreedyParams(data)
	} else if plannerName == "NewPlanner_D*" {
		params = parseNewPlannerParams(data)
	} else if plannerName == "NewPlanner_Dijkstra" {
		params = parseNewPlannerParams(data)
	}

	params["Planner"]=plannerName

	dataPoint := data[len(data)-2]
	dataPointValues := strings.Split(dataPoint, ";")

if len(dataPointValues) < 12 {
	fmt.Println(filename)
	fmt.Println(dataPoint)
}

	length := strings.TrimSpace(dataPointValues[7])
	time := strings.TrimSpace(dataPointValues[11])

	statusInt, err := strconv.ParseInt(strings.TrimSpace(dataPointValues[10]), 10, 64)
	if err != nil {
		panic(err)
	}
	if statusInt == 6 || statusInt == 5 {
		params["Solved"]="true"
	} else {
		fmt.Println(filename)
		params["Solved"]="false"
	}

	params["Solution Length"]=length
	params["Solving Time"]=time
	if precomp, ok := params["sampler_initialization_time"]; ok {
		params["Precomputation Time"] = precomp
	} else {
		params["Precomputation Time"] = "0"
	}

	return params, map[string][][]string{}, true
}