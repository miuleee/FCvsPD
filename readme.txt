Density center-based fast clustering of widefield fluorescence imaging of cortical mesoscale functional connectivity and relation to structural connectivity
=================================================
The implementation code for our recent paper:
```
Miaowen Li, Shen Gui, Qin Huang, Liang Shi, Jinling Lu, and Pengcheng Li.
"Density center-based fast clustering of widefield fluorescence imaging of cortical mesoscale functional connectivity and relation to structural connectivity". 2019.11.11. 
```

We derived the structural connection using the code (i.e., the voxel-scale model) provided by Konx et al. (https://github.com/AllenInstitute/mouse_connectivity_models). The structural connection was then compared with functional connection using Matlab code. 

============== Python Step ==================
Step 1: Download mouse_connectivity_models code (https://github.com/AllenInstitute/mouse_connectivity_models) and install :
 Python (>=2.7 or >= 3.4)
 scikit-learn (>= 0.19)
 allensdk (>= 0.14.4) 
 pynrrd (>= 0.4.0) 

Step 2: We recorded the ID number of virus tracking experiments suitable for model fitting from Allen mouse website (http://connectivity.brain-map.org/) in advance.

Step 3: Using the API (get_experiment_data) to download the data for experiments. 

Step 4: In order to fit the interpolation model, the experimental injection sites were uniformly set in the right brain (the injection sites which were originally in the left brain were flipped to the right brain accross the midline). The flipped code can be referred to the "injsite_flip.py" file

. Or use functions provided by the model API
: mcmodels.core.VoxelData (let flip_experiments=True).

Step 5: Modify configuration file. Only the structural connections of the isocortex are generated. Replace the file named "input.json" in the path "XXX\mouse_connectivity_models-master\paper". 

Step 6: Customize the experiments that participate in the fitting (see: get_voxel_data in "model_data.py"). For example,
    def get_voxel_data(self, **kwargs):
        experiment_ids2 =  [100141219,100141273,100141454,100141473,100141495,100141563,100141599,100141780,100148503,
    	100149969,112229103,112229814,112306316,112424813,112595376,112791318,112882565,112935169,
    	112951804,114290938,114292355,115956702,116903968,120875816,120916102,126852363,126861679,
    	126862385,126907302,126908007,126909424,127084296,127866392,141602484,146077302,146858755,
    	156492394,159888336,166054929,166082842,166459070,166566678,167569313,167570021,168002073,
    	168003640,168162771,175819113,176881134,177893658,181599674,181860173,182090318,182616478,
    	182803137,182814071,183618139,184167484,184169615,249327301,249396394,263781454,265289624,
    	265292679,265929196,266172624,266176167,266249483,266486371,266487079,266644610,266645328,
    	267657327,267658747,272697944,272698650,272735030,272916915,272929308,286301303,286312782,
    	286772650,287601808,288170549,288171256,292620251,293728197,294434867,294532700,294533406,
    	296052133,296052839,297233422,297592527,297652799,297654263,297668898,297669605,297670312,
    	297892130,298179622,298182842,298404154,298720191,298830161,300929973,301583889,303614706,
    	303615412,303616127,303784745,304586645,307137980,307296433,307297141,307320960,307321674,
    	307557934,307593747,307743253,309003780,309113907,479700629,479701339,479755622,479756361,
    	479980810,479981981,479982715,479983421,479984127,482578964,501115762,501484658,502140653,
    	502180994,503018656,505790715,511234957,516491813,518619451,522636747,528510546,528963991,
    	530555292,530574594,531441947,535692871,536298726,538078619,540146149,546103149,552431726,
    	557187751,559878074,560736273,562521707,566244185,566454054,577773267,583748537,584902900,
    	584903636,585025284,586042591,587659400,590548119,591612976,591622344,592540591,593018150,
    	596253012,598604288,605661910,606785720,614735393,
    	614737677,636803957,646525997,647806688,648253235,651041703,651702588,656959718,657041814]

        experiment_ids = np.unique(experiment_ids2)
        # Right brain
        data = VoxelData(self.cache, injection_structure_ids=[self.structure_id],
                         injection_hemisphere_id=2)      
        data.get_experiment_data(experiment_ids)

        return data    

Step 7: Run "run_hyperparameter_selection.py" to fit the bandwidth of the RBF kernel by employing nested cross-validation with held-out injection experiments. 

Step 8: Run "build_model.py" to derive the projection connection. (The script will generate nodes.csv.gz & weights.csv.gz).


============== Matlab Step ==================
Run "Regional_PD_VS_FC.m" to compare the functional connectivity and the structural connectivity.










