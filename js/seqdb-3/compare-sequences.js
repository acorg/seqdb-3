
function main()
{
    console.log("compare_sequences_data", compare_sequences_data);
    show_full_sequences();
}

// --------------------------------------------------------------------------------

function show_full_sequences()
{
    const tab1 = document.createElement("table");
    for (let group of compare_sequences_data.groups) {
        const tr1 = document.createElement("tr");
        tr1.innerHTML = `<td class="group-name">${group.name}</td><td><table class="group-members"></table></td>`;
        const members = tr1.querySelector(".group-members");
        for (let seq of group.seq) {
            const tr2 = document.createElement("tr");
            tr2.innerHTML = `<td class="sequence-id">${seq.id}</td><td><table class="seq-aas"><tr></tr></table></td>`;
            for (let aa of seq.seq) {
                const aa_td = document.createElement("td");
                aa_td.classList.add(`aa${aa}`);
                aa_td.innerHTML = aa;
                tr2.querySelector(".seq-aas").appendChild(aa_td)
            }
            members.appendChild(tr2);
        }
        tab1.appendChild(tr1);
    }
    document.querySelector("#full-sequences").appendChild(tab1);

    // tab1.appendChild(document.createElement("tr"
    // div.appendChild()
    // div.innerHTML = "<table></table>";
    // div.querySelector("table").appendChild("<tr><td>jopa</td></tr>");
    // div.querySelector("table").appendChild("<tr><td>jopa</td></tr>");
}

// --------------------------------------------------------------------------------

// --------------------------------------------------------------------------------

document.addEventListener("DOMContentLoaded", main);
