// NODE_PATH points to a Playwright install; URL defaults to the local preview.
const {chromium}=require('playwright');
const assert=require('node:assert/strict');
const url=process.env.CURRENTS_URL||'http://127.0.0.1:8163/';
(async()=>{
 const browser=await chromium.launch({headless:true,channel:'chrome'});
 try{
  for(const width of [1280,390]){
   const page=await browser.newPage({viewport:{width,height:900}}),errors=[];
   page.on('pageerror',e=>errors.push(String(e)));
   await page.goto(url);
   await page.waitForFunction(()=>document.querySelector('#selected').textContent.includes('illustrative current'));
   assert(await page.getByText('two alternative magnetospheric closure paths',{exact:false}).count()>0);
   for(const layer of ['Region 1 FAC sheet','R1 path 1 · direct to magnetopause','R1 path 2 · plasma sheet + boundary return','R1 · shared ionospheric closure']){
    await page.getByLabel(layer,{exact:true}).uncheck();
    await page.getByLabel(layer,{exact:true}).check();
   }
   await page.getByLabel('Tail + Chapman-Ferraro boundary return',{exact:true}).uncheck();
   await page.getByLabel('Tail + Chapman-Ferraro boundary return',{exact:true}).check();
   await page.getByLabel('BCBF · plasma flow, not current',{exact:true}).uncheck();
   await page.getByLabel('BCBF · plasma flow, not current',{exact:true}).check();
   for(const preset of ['0','1']){
    await page.locator('#preset').selectOption(preset);
    await page.getByRole('button',{name:'Global',exact:true}).click();
    const canvas=page.locator('canvas');
    const first=await canvas.screenshot();
    await page.waitForTimeout(250);
    const second=await canvas.screenshot({path:`/tmp/currents-${width}-${preset}.png`});
    assert(!first.equals(second),'animated canvas must change');
    await page.getByRole('button',{name:'Ionosphere',exact:true}).click();
    await page.getByLabel('Southern hemisphere FAC',{exact:true}).uncheck();
    await page.getByLabel('Southern hemisphere FAC',{exact:true}).check();
    const diagnostic=page.getByLabel('T96 current through GSM X–Z plane',{exact:true});
    assert.equal(await diagnostic.isChecked(),false);
    await diagnostic.check();
    assert(await page.locator('#current-colorbar').isVisible());
    await diagnostic.uncheck();
   }
   await page.getByRole('button',{name:'Pause arrows',exact:true}).click();
   await page.getByRole('button',{name:'Play arrows',exact:true}).click();
   assert.deepEqual(errors,[]);
   assert(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth));
   console.log('PASS',width,'presets, layers, hemisphere, animation, pause, layout');
   await page.close();
  }
 }finally{await browser.close();}
})().catch(e=>{console.error(e);process.exitCode=1;});
