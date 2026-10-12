const credentials = { email: 'demo@m4screen.org', password: 'M4demo2026' };
const candidates = [
  {rank:1,id:'PACERGEN01954',alias:'',drugclip:.1707,score:.9055,glide:-8.32,vina:-12.24,status:'未运行',kind:'pending'},
  {rank:2,id:'PACERGEN01409',alias:'',drugclip:.2103,score:.8913,glide:-8.54,vina:-11.37,status:'未运行',kind:'pending'},
  {rank:3,id:'PACERGEN01955',alias:'',drugclip:.1671,score:.8328,glide:-8.19,vina:-11.93,status:'未运行',kind:'pending'},
  {rank:4,id:'PACERGEN02306',alias:'',drugclip:.1904,score:.7757,glide:-8.58,vina:-10.72,status:'未运行',kind:'pending'},
  {rank:5,id:'PACERGEN00123',alias:'历史锚点 PACER0027',drugclip:.0505,score:.7307,glide:-9.58,vina:-10.51,status:'完成·待验证',kind:'complete'},
  {rank:6,id:'PACERGEN01755',alias:'历史锚点 PACER0010',drugclip:.1316,score:.7304,glide:-7.70,vina:-11.06,status:'完成·待验证',kind:'complete'},
  {rank:7,id:'PACERGEN00094',alias:'历史锚点 PACER0073',drugclip:.0763,score:.7041,glide:-8.20,vina:-10.82,status:'完成·待验证',kind:'complete'},
  {rank:8,id:'PACERGEN02353',alias:'',drugclip:.2104,score:.5936,glide:-7.18,vina:-10.48,status:'未运行',kind:'pending'}
];
const metrics = {
  roc:{base:.540,adapt:.714,gain:'+0.174',baseWidth:54,adaptWidth:71.4},
  pr:{base:.107,adapt:.155,gain:'+0.048',baseWidth:53.5,adaptWidth:77.5},
  bedroc:{base:.125,adapt:.180,gain:'+0.055',baseWidth:50,adaptWidth:72}
};

const loginPage=document.querySelector('#login-page');
const workspace=document.querySelector('#workspace');
const loginError=document.querySelector('#login-error');
function enterWorkspace(){sessionStorage.setItem('m4_demo_session','active');loginPage.classList.add('hidden');workspace.classList.add('visible');window.scrollTo(0,0)}
function leaveWorkspace(){sessionStorage.removeItem('m4_demo_session');workspace.classList.remove('visible');loginPage.classList.remove('hidden');document.querySelector('#account-menu').classList.remove('open');window.scrollTo(0,0)}
document.querySelector('#login-form').addEventListener('submit',event=>{event.preventDefault();const email=document.querySelector('#email').value.trim();const password=document.querySelector('#password').value;if(email===credentials.email&&password===credentials.password){loginError.textContent='';enterWorkspace()}else{loginError.textContent='账号或密码不正确，请使用页面中的 Demo 账号。'}});
document.querySelector('#demo-login').addEventListener('click',enterWorkspace);
document.querySelector('#logout').addEventListener('click',leaveWorkspace);
if(sessionStorage.getItem('m4_demo_session')==='active')enterWorkspace();

const accountMenu=document.querySelector('#account-menu');
document.querySelector('#account').addEventListener('click',event=>{event.stopPropagation();accountMenu.classList.toggle('open')});
document.addEventListener('click',()=>accountMenu.classList.remove('open'));
const nav=document.querySelector('.topnav nav');
document.querySelector('#mobile-nav').addEventListener('click',()=>nav.classList.toggle('open'));
document.querySelectorAll('.topnav nav a').forEach(link=>link.addEventListener('click',()=>{document.querySelectorAll('.topnav nav a').forEach(item=>item.classList.remove('active'));link.classList.add('active');nav.classList.remove('open')}));

document.querySelectorAll('.module').forEach(module=>module.querySelector('header').addEventListener('click',()=>{document.querySelectorAll('.module').forEach(item=>item.classList.remove('active'));module.classList.add('active')}));

document.querySelectorAll('[data-metric]').forEach(button=>button.addEventListener('click',()=>{const item=metrics[button.dataset.metric];document.querySelectorAll('[data-metric]').forEach(tab=>tab.classList.remove('active'));button.classList.add('active');document.querySelector('#base-value').textContent=item.base.toFixed(3);document.querySelector('#adapt-value').textContent=item.adapt.toFixed(3);document.querySelector('#gain').textContent=item.gain;document.querySelector('#base-bar').style.width=`${item.baseWidth}%`;document.querySelector('#adapt-bar').style.width=`${item.adaptWidth}%`}));

const body=document.querySelector('#candidate-body');
function renderCandidates(query=''){const normalized=query.trim().toLowerCase();const rows=candidates.filter(item=>`${item.id} ${item.alias}`.toLowerCase().includes(normalized));body.innerHTML=rows.map(item=>`<tr><td>${item.rank}</td><td>${item.id}${item.alias?`<span class="alias">${item.alias}</span>`:''}</td><td>${item.drugclip.toFixed(4)}</td><td>${item.score.toFixed(4)}</td><td>${item.glide.toFixed(2)}</td><td>${item.vina.toFixed(2)}</td><td><span class="status ${item.kind}">${item.status}</span></td></tr>`).join('')||'<tr><td colspan="7">未找到匹配候选</td></tr>'}
renderCandidates();
document.querySelector('#candidate-search').addEventListener('input',event=>renderCandidates(event.target.value));
document.querySelector('#export-csv').addEventListener('click',()=>{const header='rank,id,alias,drugclip,consensus_score,glide_bemin,vina_bemin,dynamics_status';const rows=candidates.map(item=>[item.rank,item.id,item.alias,item.drugclip,item.score,item.glide,item.vina,item.status].join(','));const blob=new Blob([[header,...rows].join('\n')],{type:'text/csv;charset=utf-8'});const url=URL.createObjectURL(blob);const link=document.createElement('a');link.href=url;link.download='m4_candidates_demo.csv';link.click();URL.revokeObjectURL(url)});

const modal=document.querySelector('#run-modal');
function openModal(){modal.classList.add('open')}
function closeModal(){modal.classList.remove('open')}
document.querySelector('#open-run').addEventListener('click',openModal);
document.querySelectorAll('[data-open-run]').forEach(button=>button.addEventListener('click',openModal));
document.querySelector('#close-run').addEventListener('click',closeModal);
modal.addEventListener('click',event=>{if(event.target===modal)closeModal()});
document.querySelector('#run-file').addEventListener('change',event=>{document.querySelector('#file-name').textContent=event.target.files[0]?.name||'支持 CSV、SDF、SMI、TXT'});
document.querySelector('#use-demo').addEventListener('click',()=>{document.querySelector('#file-name').textContent='m4_demo_2605_molecules.smi（Demo 数据）'});
document.querySelector('#run-form').addEventListener('submit',event=>{event.preventDefault();closeModal();const toast=document.querySelector('#toast');toast.classList.add('show');setTimeout(()=>toast.classList.remove('show'),2600)});
