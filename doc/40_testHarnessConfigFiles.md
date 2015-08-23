Test Harness Configuration files "run.config"
===============
These files are looked up by automated test harness systems and control the runs of a certain case.

WARNING: the user must ensure that no blank variables are specified, since the test harness might hang in such a case. The user is strongly encouraged to test the 'run.config' files before pushing them onto the test harness server.

Example content 
-----------------------------
Taken from `examples/verificationCases/NonCatalyticHeterogeneousReaction_firstStage`

````
{
  "type" : "ParScale",
  "runs" : [
    {
      "name" :          "verificationRun",
      "input_script" :  "in.parscale",
      "type" :          "ParScale/serial",
      "pre_scripts" :   ["prerun_verificationRun.sh"],
      "post_scripts" :  ["postrun_verificationRun.sh"],
      "data" : {
            "series" : [
                          {"name" : "species", "file" : "0.100000/species.json", "elements" : ["1", "2"]}
                        ],
            "plots" :   [
                          {"name" : "species", 
                          "title" : "Species", 
                          "xdata" : "VECTORLENGTH", 
                          "ydata" : ["species.1","species.2"], 
                          "xlabel" : "grid points", 
                          "ylabel" : "species", 
                          "legend" : ["Particle 1", "Particle 2"]
                          }
                        ]
      }
      
    }
  ]
}
```` 

Read next
-----------

