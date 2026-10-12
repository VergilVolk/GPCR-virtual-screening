const candidates=[
{rank:1,id:"PACERGEN01954",alias:"New candidate 01",drugclip:.1707,score:.9055,glide:-8.32,vina:-12.24,dynamic:false},
{rank:2,id:"PACERGEN01409",alias:"New candidate 02",drugclip:.2103,score:.8913,glide:-8.54,vina:-11.37,dynamic:false},
{rank:3,id:"PACERGEN01955",alias:"New candidate 03",drugclip:.1671,score:.8328,glide:-8.19,vina:-11.93,dynamic:false},
{rank:4,id:"PACERGEN02306",alias:"New candidate 04",drugclip:.1904,score:.7757,glide:-8.58,vina:-10.72,dynamic:false},
{rank:5,id:"PACERGEN00123",alias:"PACER0027",drugclip:.0505,score:.7307,glide:-9.58,vina:-10.51,dynamic:true},
{rank:6,id:"PACERGEN01755",alias:"PACER0010",drugclip:.1316,score:.7304,glide:-7.70,vina:-11.06,dynamic:true},
{rank:7,id:"PACERGEN00094",alias:"PACER0073",drugclip:.0763,score:.7041,glide:-8.20,vina:-10.82,dynamic:true},
{rank:8,id:"PACERGEN02353",alias:"New candidate 05",drugclip:.2104,score:.5936,glide:-7.18,vina:-10.48,dynamic:false}
];
const auth=document.querySelector("#auth"),workspace=document.querySelector("#workspace");
function enterWorkspace(){auth.classList.add("hidden");workspace.classList.add("active");sessionStorage.setItem("m4_demo","1")}
document.querySelector("#demo-login").addEventListener("click",enterWorkspace);
document.querySelector("#login-form").addEventListener("submit",function(e){e.preventDefault();const email=document.querySelector("#email"),password=document.querySelector("#password");if(!email.value){email.focus();return}if(!password.value){password.focus();return}enterWorkspace()});
document.querySelector("#logout").addEventListener("click",function(){sessionStorage.removeItem("m4_demo");workspace.classList.remove("active");auth.classList.remove("hidden");document.querySelector("#account-pop").classList.remove("open")});
if(sessionStorage.getItem("m4_demo")==="1")enterWorkspace();
document.querySelector("#user-menu").addEventListener("click",function(){document.querySelector("#account-pop").classList.toggle("open")});
function showPage(name){document.querySelectorAll(".app-page").forEach(function(p){p.classList.toggle("active",p.id==="page-"+name)});document.querySelectorAll(".rail nav button").forEach(function(b){b.classList.toggle("active",b.dataset.page===name)})}
document.querySelectorAll(".rail nav button").forEach(function(b){b.addEventListener("click",function(){showPage(b.dataset.page)})});
document.querySelector("[data-page-back]").addEventListener("click",function(){showPage("pipeline")});
function renderCandidates(query){query=(query||"").toLowerCase();document.querySelector("#candidate-body").innerHTML=candidates.filter(function(c){return(c.id+" "+c.alias).toLowerCase().includes(query)}).map(function(c){return'<tr data-id="'+c.id+'"><td>#'+String(c.rank).padStart(2,"0")+'</td><td><div class="candidate-name"><b>'+c.alias+'</b><small>'+c.id+'</small></div></td><td>'+c.drugclip.toFixed(4)+'</td><td><b>'+c.score.toFixed(4)+'</b></td><td>'+c.glide.toFixed(2)+'</td><td>'+c.vina.toFixed(2)+'</td><td><span class="pill '+(c.dynamic?"":"gray")+'">'+(c.dynamic?"Complete":"Not run")+'</span></td></tr>'}).join("");document.querySelectorAll("#candidate-body tr").forEach(function(row){row.addEventListener("click",function(){selectCandidate(row.dataset.id)})})}
function renderHistogram(){const min=.55,max=.95;document.querySelector("#histogram").innerHTML=candidates.map(function(c){const h=18+(c.score-min)/(max-min)*82;return'<i data-id="'+c.id+'" style="--h:'+h+'%" title="'+c.alias+' · '+c.score.toFixed(4)+'"></i>'}).join("")}
function selectCandidate(id){const c=candidates.find(function(x){return x.id===id});document.querySelectorAll("#candidate-body tr").forEach(function(r){r.classList.toggle("selected",r.dataset.id===id)});document.querySelectorAll("#histogram i").forEach(function(i){i.classList.toggle("selected",i.dataset.id===id)});document.querySelector("#selected-name").textContent=c.alias;document.querySelector("#selected-score").textContent=c.score.toFixed(4)}
document.querySelector("#candidate-search").addEventListener("input",function(e){renderCandidates(e.target.value)});
const metrics={roc:{after:"0.714",gain:"+0.174 vs baseline",base:"M48 210 C95 167 139 143 205 116 S352 72 540 18",adapt:"M48 210 C66 129 115 82 188 55 S370 28 540 18"},pr:{after:"0.155",gain:"+0.048 vs baseline",base:"M48 210 C130 194 210 174 290 145 S440 95 540 18",adapt:"M48 210 C100 174 165 148 245 118 S405 65 540 18"},bedroc:{after:"0.180",gain:"+0.055 vs baseline",base:"M48 210 C145 201 250 173 350 115 S470 45 540 18",adapt:"M48 210 C100 167 160 120 230 83 S390 34 540 18"}};
document.querySelectorAll(".metric-switch button").forEach(function(b){b.addEventListener("click",function(){document.querySelectorAll(".metric-switch button").forEach(function(x){x.classList.remove("active")});b.classList.add("active");const m=metrics[b.dataset.metric];document.querySelector("#metric-after").textContent=m.after;document.querySelector("#metric-gain").textContent=m.gain;document.querySelector(".base-line").setAttribute("d",m.base);document.querySelector(".adapted-line").setAttribute("d",m.adapt)})});
const modal=document.querySelector("#run-modal");
function openRun(){modal.classList.add("open")}
document.querySelector("#open-run").addEventListener("click",openRun);
document.querySelectorAll("[data-open-run]").forEach(function(b){b.addEventListener("click",openRun)});
document.querySelector("#close-run").addEventListener("click",function(){modal.classList.remove("open")});
modal.addEventListener("click",function(e){if(e.target===modal)modal.classList.remove("open")});
document.querySelector("#run-file").addEventListener("change",function(e){document.querySelector("#file-name").textContent=e.target.files[0]?e.target.files[0].name:"No file selected"});
document.querySelector("#use-demo").addEventListener("click",function(){document.querySelector("#file-name").textContent="m4_demo_library.csv · 2,605 structures"});
document.querySelector("#run-form").addEventListener("submit",function(e){e.preventDefault();modal.classList.remove("open");const t=document.querySelector("#toast");t.classList.add("show");setTimeout(function(){t.classList.remove("show")},3200)});
document.querySelector("#export-csv").addEventListener("click",function(){const head="rank,candidate_id,label,drugclip,consensus,glide_bemin,vina_bemin,dynamics\n",rows=candidates.map(function(c){return[c.rank,c.id,c.alias,c.drugclip,c.score,c.glide,c.vina,c.dynamic?"complete":"not_run"].join(",")}).join("\n"),a=document.createElement("a");a.href=URL.createObjectURL(new Blob([head+rows],{type:"text/csv"}));a.download="m4_candidates.csv";a.click();URL.revokeObjectURL(a.href)});
renderCandidates();renderHistogram();selectCandidate("PACERGEN01755");
